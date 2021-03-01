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
odir = f"{conf['odir']}/ggl-nside{conf['nside']}"\
       f"_{os.path.basename(conf['elledges']).replace('.txt', '')}"
print(conf, odir)

if conf['type'] == 'flask':
    real_id = int(sys.argv[2])   # Realization ID. Starts at 0
    iseed, ick = real_id // 2 + 1, real_id % 2 + 1
    print(iseed, ick)
    cshcat_fn = [f"{conf['flaskdir']}/srccat_z{iz+1}_s{iseed}_ck{ick}.parquet"
                 for iz in range(conf['nz'])]
    cshcat = [csh.cat_make_from_flaskcshcat(pd.read_parquet(fn[i]),
                                             conf['nside'])
              for fn in cshcat_fn]
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

for i in range(conf['nz_lens']):
	wc = hp.read_map(f"{conf['redmagic']}/wcountsmap_zbin{i}.fits")
	nc = hp.read_map(f"{conf['redmagic']}/ncountsmap_zbin{i}.fits")
	w = np.zeros_like(nc)
	nbar = np.zeros(len(nc))
	for z in range(0,5):
    		for i in range(len(w[z])):
        		if wc[z][i] or nc[z][i] != 0:
            			w[z][i] = wc[z][i]/nc[z][i]
		nbar[z] = sum(wc[z])/sum(w[z])
	print(f'nbar from redmagic is {nbar}')
	del(nc)

	dmap = np.zeros_like(w)
	for z in range(0,5):
    		dmap[z][w[z]>0] = wc[z][w[z]>0]/(nbar[z]*w[z][w[z]>0]) - 1
	del(wc)
	del(w)
	
	dmask = hp.read_map(f"{conf['redmagic']}/maskmap.fits")
	
	field_i = nmt.NmtField(dmask, [dmap[zi]],  purify_e=False, purify_b=False)
    for j in range(i, conf['nz_source']):
        cshcat_j = cshcat[j]
        cshmask_j = csh.mask_make(cshcat_j, conf['nside'])
        field_j = csh.field_make(cshcat_j, cshmask_j)
        w = csh.mcm_make(field_i, field_j, bins)
        cls[f'bpwrwin_{i}{j}'] = w.get_bandpower_windows()
        cls[f'cl_{i}{j}'] = w.decouple_cell(
            nmt.compute_coupled_cell(field_i, field_j))


if not os.path.exists(odir):
    os.makedirs(odir)

np.savez_compressed(ofn, **cls)
