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
    ofn = f'{odir}/cls_csh_s{iseed}_ck{ick}.npz'
elif conf['type'] == 'y1metacal':
    cshcat = mcalcat.mcalcat_process(conf['mcalcat'], conf['zbin'],
                                     conf['nside'])
    ofn = f'{odir}/cls_csh_mcal.npz'
else:
    raise ValueError(f"Computation type {conf['type']} not implemented")

bins = nmt.NmtBin.from_edges(*np.loadtxt(conf['elledges'], unpack=True,
                                         dtype='i4'))

cls = {'ell_eff': bins.get_effective_ells()}
for i in range(conf['nz']):
    cshcat_i = cshcat[i]
    cshmask_i = csh.mask_make(cshcat_i, conf['nside'])
    field_i = csh.field_make(cshcat_i, cshmask_i,
                             save_maps=conf['save_maps'],
                             maps_prefix=f'{odir}/zbin{i}')
    for j in range(i, conf['nz']):
        cshcat_j = cshcat[j]
        cshmask_j = csh.mask_make(cshcat_j, conf['nside'])
        field_j = csh.field_make(cshcat_j, cshmask_j)
        w = csh.mcm_make(field_i, field_j, bins)
        cls[f'bpwrwin_{i}{j}'] = w.get_bandpower_windows()
        cls[f'cl_{i}{j}'] = w.decouple_cell(
            nmt.compute_coupled_cell(field_i, field_j))
        if i == j:
            cls[f'nl_{i}'] = w.decouple_cell(csh.pclnoise_make(cshcat_i,
                                                               cshmask_i))
            
if not os.path.exists(odir):
    os.makedirs(odir)

print("Writing", ofn)
np.savez_compressed(ofn, **cls)
