import sys
import os
import yaml
import itertools as it
import multiprocessing as mp
import numpy as np
import healpy as hp
import pandas as pd
import pymaster as nmt
import gcl

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

# Create dirs (if needed)
if not os.path.exists(odir):
    os.makedirs(odir)

if conf['type'] == 'flask':
    real_id = int(sys.argv[2])   # Realization ID. Starts at 0
    iseed, ick = real_id // conf['nck'] + 1, real_id % conf['nck'] + 1
    print(iseed, ick)

    gclmask = f"{conf['flaskdir']}/cookies/ck{ick}.fits.gz"
    gclmask = gcl.mask_make(gclmask, ick, conf['nside'], odir)
    fsky = gclmask.mean()

    gclfield, nobj = [], []
    for iz in range(conf['nz_lns']): 
        gclcat = f"{conf['flaskdir']}/maskedcats"\
                 + f"/lnscat_z{iz+1}_s{iseed}_ck{ick}.parquet"
        gclcat = gcl.cat_fromflsk(gclcat, conf['nside'])

        f, n = gcl.field_make(gclcat, gclmask, conf['nside'],
                              save_maps=conf['save_maps'],
                              maps_prefix=f'{odir}/zbin{iz}')
        gclfield.append(f)
        nobj.append(n)
    
    ofn = f'{odir}/cls_gcl_s{iseed}_ck{ick}.npz'
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

# Pre-compute MCM WSPs (if needed)
w = gcl.mcm_make(gclmask, ick, bins, odir)

cls = {'ell_eff': bins.get_effective_ells()}
cls[f'bpwrwin'] = w.get_bandpower_windows()

# auto-correlations
for i in range(conf['nz_lns']):
    cls_coup = nmt.compute_coupled_cell(gclfield[i], gclfield[i])
    nls_coup = gcl.pclnoise_make(fsky, nobj[i], conf['nside'])
    if conf['pixwin']:
        cls_coup /= hp.pixwin(conf['nside'])**2
        nls_coup /= hp.pixwin(conf['nside'])**2
    cls[f'cl_{i}{i}'] = w.decouple_cell(cls_coup)
    cls[f'nl_{i}'] = w.decouple_cell(nls_coup)

# cross-correlations
if conf['compute_cross']:
    for i, j in it.combinations(range(conf['nz_lns']), 2):
        cls_coup = nmt.compute_coupled_cell(gclfield[i], gclfield[j])
        if conf['pixwin']:
            cls_coup /= hp.pixwin(conf['nside'])**2
        cls[f'cl_{i}{i}'] = w.decouple_cell(cls_coup)

print("Writing", ofn)
np.savez_compressed(ofn, **cls)
