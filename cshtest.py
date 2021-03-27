import sys
import os
import yaml
import itertools as it
import multiprocessing as mp
import numpy as np
import healpy as hp
import pandas as pd
import pymaster as nmt
import flask
import mcalcat
import csh

sys.stdout.flush()  # For the impatient people :).

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
    iseed, ick = real_id // conf['nck'] + 1, real_id % conf['nck'] + 1
    print(iseed, ick)
    cshcat = [f"{conf['flaskdir']}/maskedcats" +
              f"/srccat_z{iz+1}_s{iseed}_ck{ick}.parquet"
              for iz in range(conf['nz_src'])]
    cshcat = [csh.cat_fromflsk(fn, conf['nside'], conf['nonoise'])
              for fn in cshcat]
    ofn = f'{odir}/cls_csh_s{iseed}_ck{ick}.npz'
elif conf['type'] == 'y1metacal':
    cshcat = mcalcat.mcalcat_process(conf['mcalcat'], conf['zbin'],
                                     conf['nside'], conf['fgoodmap'])
    ofn = f'{odir}/cls_csh_mcal.npz'
else:
    raise ValueError(f"Computation type {conf['type']} not implemented")

# Bandpower binning - always from 0 - 3 * nside.
elledges = np.loadtxt(conf['elledges'], dtype=int)
elledges = elledges[(elledges <= 3 * conf['nside'])]
if elledges[0] > 0:
    elledges = np.insert(elledges, 0, 0)
if elledges[-1] < 3 * conf['nside']:
    elledges = np.append(elledges, 3 * conf['nside'])
bins = nmt.NmtBin.from_edges(elledges[:-1], elledges[1:])

cls = {'ell_eff': bins.get_effective_ells()}
for i in range(conf['nz_src']):
    cshcat_i = cshcat[i]
    cshmask_i = csh.mask_make(cshcat_i, conf['nside'])
    field_i = csh.field_make(cshcat_i, cshmask_i,
                             save_maps=conf['save_maps'],
                             maps_prefix=f'{odir}/zbin{i}')
    for j in range(i, conf['nz_src']):
        cshcat_j = cshcat[j]
        cshmask_j = csh.mask_make(cshcat_j, conf['nside'])
        field_j = csh.field_make(cshcat_j, cshmask_j)
        w = csh.mcm_make(field_i, field_j, bins)
        cls[f'bpwrwin_{i}{j}'] = w.get_bandpower_windows()
        cls_coup = nmt.compute_coupled_cell(field_i, field_j)
        if conf['pixwin']:
            cls_coup /= np.array([hp.pixwin(conf['nside'])] * 4)**2
        cls[f'cl_{i}{j}'] = w.decouple_cell(cls_coup)
        if i == j:
            nls_coup = csh.pclnoise_make(cshcat_i, cshmask_i)
            if conf['pixwin']:
                nls_coup /= np.array([hp.pixwin(conf['nside'])] * 4)**2
            cls[f'nl_{i}'] = w.decouple_cell(nls_coup)
            
if not os.path.exists(odir):
    os.makedirs(odir)

print("Writing", ofn)
np.savez_compressed(ofn, **cls)
