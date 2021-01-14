#     Contamined shear maps!
#       cookie and seed independence  
#       

import numpy as np
import healpy as hp
import pymaster as nmt
import time
import getpass;
import sys
import os

USER=getpass.getuser()

# Avoid to measure twice
overwrite_cls=False
try:
    sd= sys.argv[1]# VARIABLE_FROM_NERSC
    ick = sys.argv[2]
    print( f"Starting  seed: {sd} cookie: {ick}...")
except Exception as err:
    print("Error in argv: ", err)
    exit()


root = '/global/cscratch1/sd/faoli/data-des-y1/'
outcls   = '/global/cscratch1/sd/faoli/outggl/'



fracdet_min = 0.0
NSIDE = 2048
pur_e = False
pur_b = False

if USER == 'foliveira':
    dryrun  = True
    NSIDE=128
    root =  '/home/foliveira/programs/data-des-y1/seed_hugo/dryrun/'
    outcls = f'{root}/outggl/'

dmapfn   = lambda sd, zbin, NSIDE: f'{root}map/dmap_z{zbin}_nside{NSIDE}_s{sd}.fits'
ggmapfn  = lambda sd, zbin, NSIDE: f'{root}sn/gg-s{sd}-f2z{zbin}_nside{NSIDE}_sn.fits.gz'
maskfn   = lambda ck, NSIDE: f'{root}mask/ck{ck}_nside{NSIDE}.fits'


print("Loading binning scheme...")
binnfn = f'{root}binning/bin_log_N20_lmin30_lmax3000.dat'
bpws, ells, weights = np.loadtxt(binnfn, unpack=True)
b = nmt.NmtBin(
    NSIDE,
    bpws=bpws.astype('i4'),
    ells=ells.astype('i4'),
    weights=weights.astype('f4'))
elllog = b.get_effective_ells()





for ck in [ick]: #range(1,9):
    print(f"cookie {ck}")
    mask = hp.read_map(maskfn(ck, NSIDE), dtype=['f4'], verbose=False)
    mask[(mask <= fracdet_min)] = 0.0
    print(f"= f_sky: {mask.sum()/hp.nside2npix(NSIDE)}")

    print("Creating working space...")
    w = nmt.NmtWorkspace()
    w.read_from(f"{root}coupling/gglNoPur-ck{ck}_nside{NSIDE}_bin_log_N20_lmin30_lmax3000.coup")

    t0 =time.time()
    for zj in range(1,5):
        t1 =time.time()

        ggmap = hp.read_map(ggmapfn(sd, zj, NSIDE), field=[0,1], dtype=['f4']*2, verbose=False)

        print("Nmt Field gamma")
        fl2 = nmt.NmtField(mask, [ggmap[0], ggmap[1]],
                    purify_e=pur_e, purify_b=pur_b,
                    n_iter=3 ) #, n_iter_mask_purify=10)

        for zi in range(1,6):
            if not overwrite_cls:
                if os.path.isfile(f"{outcls}cl_sn_NoPur-s{sd}_ck{ck}_f1z{zi}f2z{zj}_nside{NSIDE}_bin_log_N20_lmin30_lmax3000.dat"):
                    print(f"skiping s{sd}_c{ck} f1z{zi} f2z{zj}")
                    break

            dmap = hp.read_map(dmapfn(sd, zi, NSIDE), dtype=['f4'], verbose=False)

            print(f"Preparing field f1z{zi} f2z{zj}...")
            fl0 = nmt.NmtField(mask, [dmap],
                purify_e=pur_e, purify_b=pur_b,
                n_iter=3) #, n_iter_mask_purify=10)

            cl_nmt = nmt.compute_full_master(fl0, fl2, b, workspace=w)
            print("saving cls...")
            np.savetxt(f"{outcls}cl_sn_NoPur-s{sd}_ck{ck}_f1z{zi}f2z{zj}_nside{NSIDE}_bin_log_N20_lmin30_lmax3000.dat",
                         np.c_[elllog, cl_nmt[0], cl_nmt[1]], header="ell \t gal-E \t gal-B")

            if False: # sd==1 and ck == 1:
                print("Cl Full Sky (for sampling)")
                cl_fs = hp.anafast([dmap, ggmap[0], ggmap[1]])
                np.savetxt(f"{outcls}cl_FS-s{sd}_ck{ck}_f1z{zi}f2z{zj}_nside{NSIDE}_bin_log_N20_lmin30_lmax3000.dat",
                            cl_fs, header="# (TT, EE, BB, TE, EB, TB)")
        print(f"Looping time {time.time() - t1:3.2f}s")
print(f'Execution time (cumulative) f2z{zj} f1z1-5): {time.time()-t0: 3.2}s' )
