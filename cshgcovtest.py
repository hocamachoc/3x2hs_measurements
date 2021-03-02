import sys
import os
import yaml
import itertools as it
import multiprocessing as mp
import numpy as np
import pandas as pd
import pymaster as nmt
import flask
import mcalcat
import csh


# Configuration file
conf = sys.argv[1]
with open(conf, 'rt') as f:
    conf = yaml.safe_load(f)
odir = f"{conf['odir']}/nside{conf['nside']}"\
       f"_{os.path.basename(conf['elledges']).replace('.txt', '')}"
if conf['nonoise']:
    odir += '_nonoise'
print(conf, odir)

if conf['type'] == 'flask':
    real_id = int(sys.argv[2])   # Realization ID. Starts at 0
    iseed, ick = real_id // 2 + 1, real_id % 2 + 1
    print(iseed, ick)
    cshcat = [f"{conf['flaskdir']}/srccat_z{iz+1}_s{iseed}_ck{ick}.parquet"
              for iz in range(conf['nz'])]
    cshcat = [csh.cat_fromflsk(fn, conf['nside'], conf['nonoise'])
              for fn in cshcat]
    print(cshcat)
    ofn = f'{odir}/gcov_csh_s{iseed}_ck{ick}.npz'
elif conf['type'] == 'y1metacal':
    cshcat = mcalcat.mcalcat_process(conf['mcalcat'], conf['zbin'],
                                     conf['nside'])
    ofn = f'{odir}/gcov_csh_mcal.npz'
else:
    raise ValueError(f"Computation type {conf['type']} not implemented")

bins = nmt.NmtBin.from_edges(*np.loadtxt(conf['elledges'], unpack=True,
                                         dtype='i4'))

print("Constructing fields")
cshmask = [csh.mask_make(cshcat[i], conf['nside'])
           for i in  range(conf['nz'])]
field = [csh.field_make(cshcat[i], cshmask[i])
         for i in range(conf['nz'])]

pairs = [(i, j)
         for i, j in it.combinations_with_replacement(range(conf['nz']), 2)]

print("Preparing PCL workspaces")
w = [csh.mcm_make(field[i], field[j], bins) for i, j in pairs]

print("Loading input Cls")
lmax = 3 * conf['nside']
cl_in = [[None for i in range(conf['nz'])] for j in range(conf['nz'])]
for ii, (i, j) in enumerate(pairs):
    cl = flask.load_inputcl(i, j, lmax)
    print(cl.shape, w[ii].wsp.lmax+1, type(w[ii].wsp.lmax+1), ii, type(ii))
    cl[:, :w[ii].wsp.lmax+1] = w[ii].couple_cell(cl)
    cl[:, w[ii].wsp.lmax:] = cl[:, w[ii].wsp.lmax][:, None]
    if i == j:
        cl += csh.pclnoise_make(cshcat[i], cshmask[i])
    cl /= np.mean(cshmask[i] * cshmask[j])
    cl_in[i][j] = cl
    if i != j:
        cl_in[j][i] = cl
print(len(cl_in), len(cl_in[0]), len(cl_in[0][0]), cl_in[0][0][0].shape)

print("Building GCOV")
n_ell = bins.get_n_bands()
cov = {}
for a in range(len(pairs)):
    for b in range(a, len(pairs)):
        a1, a2 = pairs[a]
        b1, b2 = pairs[b]
        print(f'{a1}{a2}_{b1}{b2}')
        cov[f'{a1}{a2}_{b1}{b2}'] = csh.gcov_make(field[a1], field[a2],
                                                  field[b1], field[b2],
                                                  w[a], w[b],
                                                  cl_in[a1][b1], cl_in[a1][b2],
                                                  cl_in[a2][b1], cl_in[a2][b2],
                                                  n_ell)
            
if not os.path.exists(odir):
    os.makedirs(odir)

print("Writing", ofn)
np.savez_compressed(ofn, **cov)
