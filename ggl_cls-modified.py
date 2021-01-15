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

# USER=getpass.getuser()

# Avoid to measure twice
# # overwrite_cls=False
# # try:
# #     sd= sys.argv[1]# VARIABLE_FROM_NERSC
# #     ick = sys.argv[2]
# #     print( f"Starting  seed: {sd} cookie: {ick}...")
# # except Exception as err:
# #     print("Error in argv: ", err)
# #     exit()


# root = '/global/cscratch1/sd/faoli/data-des-y1/'
# outcls   = '/global/cscratch1/sd/faoli/outggl/'



fracdet_min = 0.0
NSIDE = 1024
pur_e = False
pur_b = False

# if USER == 'foliveira':
#     dryrun  = True
#     NSIDE=128
#     root =  '/home/foliveira/programs/data-des-y1/seed_hugo/dryrun/'
#     outcls = f'{root}/outggl/'

# dmapfn   = lambda sd, zbin, NSIDE: f'{root}map/dmap_z{zbin}_nside{NSIDE}_s{sd}.fits'
# ggmapfn  = lambda sd, zbin, NSIDE: f'{root}sn/gg-s{sd}-f2z{zbin}_nside{NSIDE}_sn.fits.gz'
# maskfn   = lambda ck, NSIDE: f'{root}mask/ck{ck}_nside{NSIDE}.fits'

maskfn  = lambda zbin: f'maps/metacal_maps_nside1024/maskmap_zbin{zbin}.fits'
wmapfn  = lambda zbin: f'maps/redmagic_maps_nside1024/wcountsmap_zbin{zbin}.fits'
nmapfn  = lambda zbin: f'maps/redmagic_maps_nside1024/countsmap_zbin{zbin}.fits'
ggmapfn = lambda i, zbin: f'maps/metacal_maps_nside1024/g{i}map_zbin{zbin}.fits'


print("Loading binning scheme...")
binnfn = f'datay1/bin_log_N20_lmin30_lmax3000.dat'
bpws, ells, weights = np.loadtxt(binnfn, unpack=True)
b = nmt.NmtBin(
    NSIDE,
    bpws=bpws.astype('i4'),
    ells=ells.astype('i4'),
    weights=weights.astype('f4'))
elllog = b.get_effective_ells()





# for ck in [ick]: #range(1,9):
#     print(f"cookie {ck}")
#     mask = hp.read_map(maskfn(ck, NSIDE), dtype=['f4'], verbose=False)
#     mask[(mask <= fracdet_min)] = 0.0
#     print(f"= f_sky: {mask.sum()/hp.nside2npix(NSIDE)}")

#     print("Creating working space...")
#     w = nmt.NmtWorkspace()
#     w.read_from(f"{root}coupling/gglNoPur-ck{ck}_nside{NSIDE}_bin_log_N20_lmin30_lmax3000.coup")

# t0 =time.time()
# for zj in range(1,5):
#     t1 =time.time()
    
mask = hp.read_map(maskfn(1), dtype=['f4'], verbose=False)

g1map = hp.read_map(ggmapfn(1,1), field=[0], dtype=['f4'], verbose=False)
g2map = hp.read_map(ggmapfn(2,1), field=[0], dtype=['f4'], verbose=False)
ggmap = [g1map,g2map]
ggmap = np.array(ggmap)
print(ggmap.shape)

print("Nmt Field gamma")
fl2 = nmt.NmtField(mask, [ggmap[0], ggmap[1]],
            purify_e=pur_e, purify_b=pur_b,
            n_iter=3 ) #, n_iter_mask_purify=10)

#     for zi in range(1,6):
#         if not overwrite_cls:
#             if os.path.isfile(f"{outcls}cl_sn_NoPur-s{sd}_ck{ck}_f1z{zi}f2z{zj}_nside{NSIDE}_bin_log_N20_lmin30_lmax3000.dat"):
#                 print(f"skiping s{sd}_c{ck} f1z{zi} f2z{zj}")
#                 break

wc = hp.read_map(wmapfn(1), dtype=['f4'], verbose=False)
nc = hp.read_map(nmapfn(1), dtype=['f4'], verbose=False)

w = np.zeros(len(nc))
for i in range(len(w)):
    if wc[i] or nc[i] != 0:
        w[i] = wc[i]/nc[i]
print(len(w))
nbar = sum(wc)/sum(w)

dmap = np.zeros(len(wc))
dmap[w>0] = wc[w>0]/(nbar*w[w>0]) - 1

# print(f"Preparing field f1z{zi} f2z{zj}...")
fl0 = nmt.NmtField(mask, [dmap],
    purify_e=pur_e, purify_b=pur_b,
    n_iter=3) #, n_iter_mask_purify=10)

cl_nmt = nmt.compute_full_master(fl0, fl2, b) #, workspace=w)

#         print("saving cls...")
#         np.savetxt(f"{outcls}cl_sn_NoPur-s{sd}_ck{ck}_f1z{zi}f2z{zj}_nside{NSIDE}_bin_log_N20_lmin30_lmax3000.dat",
#                      np.c_[elllog, cl_nmt[0], cl_nmt[1]], header="ell \t gal-E \t gal-B")

#         if False: # sd==1 and ck == 1:
#             print("Cl Full Sky (for sampling)")
#             cl_fs = hp.anafast([dmap, ggmap[0], ggmap[1]])
#             np.savetxt(f"{outcls}cl_FS-s{sd}_ck{ck}_f1z{zi}f2z{zj}_nside{NSIDE}_bin_log_N20_lmin30_lmax3000.dat",
#                         cl_fs, header="# (TT, EE, BB, TE, EB, TB)")
#     print(f"Looping time {time.time() - t1:3.2f}s")
# print(f'Execution time (cumulative) f2z{zj} f1z1-5): {time.time()-t0: 3.2}s' )



##############################################################################
############# UNCOMMENT THIS SECTION IF YOU WANT TO PLOT #####################

#import matplotlib.pyplot as plt

#plt.loglog(elllog, cl_nmt[0], label='TE', marker='o')
#plt.loglog(elllog, cl_nmt[1], label='TB', marker='o')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('$\\ell$', fontsize=16)
#plt.ylabel('$C_\\ell$', fontsize=16)
#plt.legend(loc='upper right', ncol=2, labelspacing=0.1)
#plt.title('galaxy-galaxy lensing')
#plt.show()
