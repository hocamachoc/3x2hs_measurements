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
import csh, gcl
from functools import partial
from os.path import exists, getsize

print = partial(print, flush=True)  # For the impatient people :).

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

def mcm_process(ick, zj, dmask_i, cshmask_j, bins, mcm_dir):
    """ Return MCM workspace for pos - pos
    """
    print(f'calculating mcm - ck{ick} ...', flush=True)
    path = f"{mcm_dir}/mcm_ggl_ck{ick}_source_z{zj}.fits"
    w = nmt.NmtWorkspace()
    if exists(path) and getsize(path) > 0:
        w.read_from(path)
        print(f'mcm file already exists', flush=True)
        return w
    print(f"mcm file doesn't exist", flush=True)
    f0 = nmt.NmtField(dmask_i, [dmask_i])
    f2 = nmt.NmtField(cshmask_j, [cshmask_j, cshmask_j])
    w.compute_coupling_matrix(f0, f2, bins)
    w.write_to(path)
    print(f'mcm file computed and saved at {path}', flush=True)
    return w

if conf['type'] == 'flask':
    real_id = int(sys.argv[2])   # Realization ID. Starts at 0
    iseed, ick = real_id // conf['nck'] + 1, real_id % conf['nck'] + 1
    print(iseed, ick)

    # Prepare cosmic-shear stuff
    cshcat = [f"{conf['flaskdir']}/maskedcats" +
              f"/srccat_z{iz+1}_s{iseed}_ck{ick}.parquet"
              for iz in range(conf['nz_src'])]
    cshcat = [csh.cat_fromflsk(fn, conf['nside'], conf['nonoise'])
              for fn in cshcat]

    # Prepare galaxy-clustering stuff
    gclmask = f"{conf['flaskdir']}/cookies/ck{ick}.fits.gz"
    gclmask = gcl.mask_make(gclmask, ick, conf['nside'], odir)
    fsky = gclmask.mean()
    gclfield, nobj = [], []
    print(f'creating clustering field',flush=True)
    for iz in range(conf['nz_lns']): 
        print(f'zbin {iz}',flush=True)
        gclcat = f"{conf['flaskdir']}/maskedcats"\
                 + f"/lnscat_z{iz+1}_s{iseed}_ck{ick}.parquet"
        gclcat = gcl.cat_fromflsk(gclcat, conf['nside'])

        f, n = gcl.field_make(gclcat, gclmask, conf['nside'],
                              save_maps=conf['save_maps'],
                              maps_prefix=f'{odir}/zbin{iz}')
        gclfield.append(f)
        nobj.append(n)

    # Output filename
    ofn = f'{odir}/cls_ggl_s{iseed}_ck{ick}.npz'
else:
    raise NotImplementedError(f"Computation type {conf['type']} not"
                              + " implemented")

# Bandpower binning - always from 0 to 3 * nside.
elledges = np.loadtxt(conf['elledges'], dtype=int)
elledges = elledges[(elledges <= 3 * conf['nside'])]
if elledges[0] > 0:
    elledges = np.insert(elledges, 0, 0)
if elledges[-1] < 3 * conf['nside']:
    elledges = np.append(elledges, 3 * conf['nside'])
bins = nmt.NmtBin.from_edges(elledges[:-1], elledges[1:])

cls = {'ell_eff': bins.get_effective_ells()}

if conf['nside'] <= 2048:
    field_i = gclfield
    field_j = []
    print(f'creating shear field',flush=True)
    for j in range(conf['nz_src']):
        print(f'zbin {j}',flush=True)
        cshcat_j = cshcat[j]
        cshmask_j = csh.mask_make(cshcat_j, conf['nside'])
        field_j.append(csh.field_make(cshcat_j, cshmask_j,
                                 save_maps=conf['save_maps'],
                                 maps_prefix=f'{odir}/zbin{j}'))

    for i in range(conf['nz_lns']):
        for j in range(conf['nz_src']):
            print(f'computing cls z{i+1} z{j+1}', flush=True)
            #w.compute_coupling_matrix(field_i[i], field_j[j], bins)
            w = mcm_process(conf['nck'], j+1, gclmask, cshmask_j[j], bins, conf['mcm_dir'])
            cls[f'bpwrwin_{i}{j}'] = w.get_bandpower_windows()
            cls_coup = nmt.compute_coupled_cell(field_i[i], field_j[j])
            if conf['pixwin']:
                cls_coup /= np.array([hp.pixwin(conf['nside'])] * 2)**2
            cls[f'cl_{i}{j}'] = w.decouple_cell(cls_coup)
            print('\n', flush=True)

else: 
    for i in range(conf['nz_lns']):
        field_i = gclfield[i]
        for j in range(conf['nz_src']):
            cshcat_j = cshcat[j]
            cshmask_j = csh.mask_make(cshcat_j, conf['nside'])
            field_j = csh.field_make(cshcat_j, cshmask_j,
                                 save_maps=conf['save_maps'],
                                 maps_prefix=f'{odir}/zbin{j}')
            w.compute_coupling_matrix(field_i, field_j, bins)
            cls[f'bpwrwin_{i}{j}'] = w.get_bandpower_windows()
            cls_coup = nmt.compute_coupled_cell(field_i, field_j)
            if conf['pixwin']:
                cls_coup /= np.array([hp.pixwin(conf['nside'])] * 2)**2
            cls[f'cl_{i}{j}'] = w.decouple_cell(cls_coup)

print("Writing", ofn)
np.savez_compressed(ofn, **cls)
