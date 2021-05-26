import sys
import os
import yaml
import itertools as it
import multiprocessing as mp
import numpy as np
import pandas as pd
import pymaster as nmt
import healpy as hp
import flask
import mcalcat
import csh
from os.path import exists, getsize


# Configuration file
conf = sys.argv[1]
with open(conf, 'rt') as f:
    conf = yaml.safe_load(f)
odir = f"{conf['odir']}/ggl-nside{conf['nside']}"\
       f"_{os.path.basename(conf['elledges']).replace('.txt', '')}"
print(conf, odir)

bins = nmt.NmtBin.from_edges(*np.loadtxt(conf['elledges'], unpack=True,
                                         dtype='i4'))

cls = {'ell_eff': bins.get_effective_ells()}

def mcm_process(ick, zj, dmask_i, cshmask_j, bins, mcm_dir):
    """ Return MCM workspace for pos - pos
    """
    print(f'mcm creation - ck{ick} source{zj}: start', flush=True)
    path = f"{mcm_dir}/mcm_ggl_ck{ick}_source_z{zj}.fits"
    w = nmt.NmtWorkspace()
    if exists(path) and getsize(path) > 0:
        w.read_from(path)
        print(f'mcm creation - ck{ick} source{zj}: already exists',flush=True)
        return w
    f0 = nmt.NmtField(dmask_i, [dmask_i])
    f2 = nmt.NmtField(cshmask_j, [cshmask_j, cshmask_j])
    w.compute_coupling_matrix(f0, f2, bins)
    w.write_to(path)
    print('mcm creation - ck{ick} source{zj}: created',flush=True)
    return w

if conf['type'] == 'flask':
    #real_id = int(sys.argv[2])   # Realization ID. Starts at 0
    #iseed, ick = real_id // 2 + 1, real_id % 2 + 1
    #print(iseed, ick)
    
    iseed = 1
    ick = int(sys.argv[2])
    
    #    print(f"\n SEED {iseed} COOKIE {ick} BEGAN \n")
    print(f"processing cshcat ...")
    cshcat_fn = [f"{conf['flaskdir']}/srccat_z{iz+1}_s{iseed}_ck{ick}.parquet"
                     for iz in range(conf['nz_source'])]
    cshcat = [csh.cat_fromflsk(fn, conf['nside'])
                  for fn in cshcat_fn]
    #    ofn = f'{odir}/cls_ggl_s{iseed}_ck{ick}.npz'
    
    #    field_i = []
    print(f"getting dmask ...")
    dmask_i = hp.read_map(f'{conf["ck_dir"]}/ck{ick}_nside{conf["nside"]}.fits')
    #    print(f"getting density fields ...")
    #    for i in range(conf['nz_lens']):
    #        dmap_i = hp.read_map(f'{conf["dmap_dir"]}/dmap_z{i}_nside{conf["nside"]}_s{iseed}.fits')	
    #        field_i.append(nmt.NmtField(dmask_i, [dmap],  purify_e=False, purify_b=False))

    #    field_j = []
    cshmask_j = []
    print(f"getting shear mask")
    for j in range(conf['nz_source']):
        cshcat_j = cshcat[j]
        cshmask_j.append(csh.mask_make(cshcat_j, conf['nside']))
    #        field_j.append(csh.field_make(cshcat_j, cshmask_j))
    
    #for i in range(conf['nz_lens']):
    for j in range(conf['nz_source']):
        w = mcm_process(ick, j+1, dmask_i, cshmask_j[j], bins, conf['mcm_dir'])
            #    cls[f'bpwrwin_{i}{j}'] = w.get_bandpower_windows()
            #    cls[f'cl_{i}{j}'] = w.decouple_cell(
            #        nmt.compute_coupled_cell(field_i, field_j))


elif conf['type'] == 'y1metacal':
    cshcat = mcalcat.mcalcat_process(conf['mcalcat'], conf['zbin'],
                                     conf['nside'])
    ofn = f'{odir}/cls_ggl_mcal.npz'
    
    for i in range(conf['nz_lens']):
        wc = hp.read_map(f"{conf['redmagic']}/wcountsmap_zbin{i}.fits")	
        dmask = hp.read_map(f"{conf['redmagic']}/maskmap.fits")
    
        nbar = sum(wc[dmask>0])/sum(dmask[dmask>0])
        print(nbar)

        dmap = np.full(len(dmask), 0.0)
        dmap[dmask>0] = wc[dmask>0]/(nbar*dmask[dmask>0]) - 1
	
        field_i = nmt.NmtField(dmask, [dmap],  purify_e=False, purify_b=False)
        for j in range(conf['nz_source']):
            print(i,j)
            cshcat_j = cshcat[j]
            cshmask_j = csh.mask_make(cshcat_j, conf['nside'])
            field_j = csh.field_make(cshcat_j, cshmask_j)
            w = csh.mcm_make(field_i, field_j, bins)
            cls[f'bpwrwin_{i}{j}'] = w.get_bandpower_windows()
            cls[f'cl_{i}{j}'] = w.decouple_cell(
                nmt.compute_coupled_cell(field_i, field_j))
else:
    raise ValueError(f"Computation type {conf['type']} not implemented")


#if not os.path.exists(odir):
#    os.makedirs(odir)

#print(f"saving cls ...")
#np.savez_compressed(ofn, **cls)

#print('DONE \n')
