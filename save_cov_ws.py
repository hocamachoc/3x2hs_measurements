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
import gcl
import time
import multiprocessing as mp

t0 = time.time()

sys.stdout.flush()  # For the impatient people :).

#-----------------------------------------------
# Defining some functions 

def save_cov_workspace(fa1, fa2, fb1, fb2, fn):
    if not os.path.exists(fn):
        print(f'There is no {fn}\n')
        print(f'Creating covariance workspace')

        cw = nmt.NmtCovarianceWorkspace()
        cw.compute_coupling_coefficients(fa1,
                                         fb1,
                                         fa2,
                                         fb2,)
        cw.write_to(fn)
        return print(f'Covariance workspace saved at {fn}')

    else:
        return print(f'Covariance workspace ALREADY EXISTS {fn}')

#-----------------------------------------------
# Defining paths

data_path = "/global/cscratch1/sd/ljfaga/DESY3_data/"
theo_path = "/global/homes/l/ljfaga/3x2hs_measurements/conf/Cl_flaskv2p0_nolimber_emu_Nsource4/"

print(f'Catalogs path: {data_path}')
print(f'Theoretical C_ells path: {theo_path}\n')

#-----------------------------------------------
# Configuration file
print('Reading configuration file \n')

conf = sys.argv[1]
#conf = 'etc/y3data-LJF.yml'
with open(conf, "rt") as f:
    conf = yaml.safe_load(f)
odir = (
    f"{conf['odir']}/nside{conf['nside']}"
    f"_{os.path.basename(conf['elledges']).replace('.txt', '')}"
)

print(f'Printing configuration file info:\n {conf} \n')

block = sys.argv[2]

print(f'Chosen covariance block: {block}')

if block != 'ggl-ggl' and block != 'gcl-ggl' and block != 'gcl-gcl': 
    raise Exception(f"Couldn't understand {block}." + 
                "The second argument of this script should be one of the following:" +
                 "'ggl-ggl', 'gcl-ggl', or 'gcl-gcl'")

comb = int(sys.argv[3])

print(f'Calculating and saving covariance workspaces for block {block} combination #{comb}\n')

#----------------------------------------------
# Getting bandpower binning scheme (always from 0 - 3 * nside)
print('Getting bandpower binning scheme \n')

elledges = np.loadtxt(conf["elledges"], dtype=int)
elledges = elledges[(elledges <= 3 * conf["nside"])]
if elledges[0] > 0:
    elledges = np.insert(elledges, 0, 0)
if elledges[-1] < 3 * conf["nside"]:
    elledges = np.append(elledges, 3 * conf["nside"])
bins = nmt.NmtBin.from_edges(elledges[:-1], elledges[1:])

#----------------------------------------------
# Getting shape masks and fields                                                                                                
print('Getting galaxy shapes masks and fields\n')

if conf["type"] == "y3data":
    # Prepare cosmic-shear stuff
    cshcat_full = csh.cat_fromy3data(conf["metacal"], conf["nside"])
    cshcat = [
        cshcat_full.loc[cshcat_full["bin_number"] == iz]
        for iz in range(conf["nz_src"])
    ]

cshmask = [
    csh.mask_make(cshcat[i], conf["nside"]) for i in range(conf["nz_src"])
]
cshfield = [csh.field_make(cshcat[i], cshmask[i]) for i in range(conf["nz_src"])]

#---------------------------------------------
# Getting galaxy position mask and fields
print('Getting galaxy position mask and fields\n')

if conf["type"] == "y3data":
        # Prepare galaxy-clustering stuff
    gclmask = gcl.mask_make_y3data(conf["redmagic_mask"], conf["nside"])
    fsky = gclmask.mean()
    gclfield, nobj = [], []
    gclcat_full = gcl.cat_fromy3data(conf["redmagic"], conf["nside"])
    for iz in range(conf["nz_lns"]):
        gclcat = gclcat_full.loc[gclcat_full["bin_number"] == iz + 1]

        f, n = gcl.field_make(
            gclcat,
            gclmask,
            conf["nside"],
            save_maps=conf["save_maps"],
            maps_prefix=f"{odir}/zbin{iz}",
        )
        gclfield.append(f)
        nobj.append(n)

#--------------------------------------------
# Defining zbins combinations for csh, gcl, and ggl

csh_pairs = [
    (i, j)
    for i, j in it.combinations_with_replacement(range(conf["nz_src"]), 2)
]

gcl_pairs = [
    (i, i)
    for i in range(5)
]

ggl_pairs = [
    (i, j)
    for i in range(5) for j in range(4)
]

print(f'Preparing fields took {(time.time() - t0)/60} minutes')

print(f'Saving covariance workspaces for block {block}\n')

cws_path = '/global/cscratch1/sd/ljfaga/DESY3_data/cov_workspace/'
print(f'Saving them at {cws_path} \n')


#---#---#

if block == 'gcl-ggl':

    print('Clustering-GGL block')

    covpairs = [(gcl_pairs[a], ggl_pairs[b])
             for a in range(len(gcl_pairs))
             for b in range(a, len(ggl_pairs))]

    covpairs = covpairs[comb][:]

    a1, a2 = covpairs[0]
    b1, b2 = covpairs[1]

    fn = f'{cws_path}/gcl-ggl_cws_{a1}{a2}_{b1}{b2}'

    t0 = time.time()

    save_cov_workspace( gclfield[a1],
                        gclfield[a2],
                        gclfield[b1],
                        cshfield[b2],
                        fn)


    print(f'total time = {(time.time() - t0 )/60} minutes')

#---#---#

if block == 'ggl-ggl':

    print('GGL-GGL block')

    covpairs = [(ggl_pairs[a], ggl_pairs[b])
             for a in range(len(ggl_pairs))
             for b in range(a, len(ggl_pairs))]
    
    covpairs = covpairs[comb][:]

    a1, a2 = covpairs[0]
    b1, b2 = covpairs[1]

    fn = f'{cws_path}/ggl-ggl_cws_{a1}{a2}_{b1}{b2}'

    t0 = time.time()

    save_cov_workspace( gclfield[a1],
                        gclfield[b1],
                        cshfield[a2],
                        cshfield[b2],
                        fn)


    print(f'total time = {(time.time() - t0 )/60} minutes')

print('Finished!!!')
