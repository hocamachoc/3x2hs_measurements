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

sys.stdout.flush()  # For the impatient people :).

# TO DO
# GAUSSIAN CODE:
# - perguntar para o hugo como ele fez a csh gcov 
# (usou flasko ou y1data para o field? ou o código não está atualizado? )
# - calculate or import every workspace
# - add clustering noise
# - add loop for all zbins combinations
# - figure out how i should calculate cross csh, ggl, gcl terms
# - validate with previous covariance
# - compute cov coef only once? (as done in namaster cov tutorial)
# - so far code only works for y3data

#-----------------------------------------------
# Defining some functions 
# TO DO: organize it better (like what was done for shear)

# Make clustering coupling matrix workspace
def mcm_make(mask_fn, mask, b, cachedir):
    """
    Returns MCM workspace for position-position (clustering)
    """
    path = mask_fn.replace(".fits", "_wspmcm.fits")
    path = path.replace(".gz", "")
    w = nmt.NmtWorkspace()
    #if os.path.exists(path) and os.path.getsize(path) > 0:
    #    w.read_from(path)
    #    return w
    f1 = nmt.NmtField(mask, [mask])
    w.compute_coupling_matrix(f1, f1, b)
    w.write_to(path)
    return w
    
# Read an existing workspace
def mcm_read(path):
    w = nmt.NmtWorkspace()
    w.read_from(path)
    return w
    
# Read theoretical cls
def load_inputcl(i, j, lmax, idir, cl_type):
    """
    Loads FLASK input Cls to a numpy array
    """
    if cl_type=='csh':
        fn = f'{idir}/Y3_5x2pt_Nsource4-Cl_f10z{i}f10z{j}.dat'
    if cl_type=='gcl':
        fn = f'{idir}/Y3_5x2pt_Nsource4-Cl_f{i}z{i}f{j}z{j}.dat'
    if cl_type=='ggl':
        fn = f'{idir}/Y3_5x2pt_Nsource4-Cl_f10z{j}f{i}z{i}.dat'
    # TODO: Introduce a check here. Normally, l starts at 1 on FLASK/CLike
    # predictions, but we can introduce a check to allow for generality.
    cl = np.loadtxt(fn, usecols=1)
    cl = np.insert(cl, 0, 0)[:lmax]
    
    if cl_type=='csh':
        return np.array([cl] + [np.zeros_like(cl)] * 3)
    if cl_type=='gcl':
        return np.array([cl])
    if cl_type=='ggl':
        return np.array([cl] + [np.zeros_like(cl)] )

# Calculate clustering-clustering gaussian covariance
def gcl_gcov_make(
    fa1, fa2, fb1, fb2, wa, wb, cla1b1, cla1b2, cla2b1, cla2b2, n_ell, cws_path
):
    """Returns the Gaussian covariance matrix fa1fa2_fb1fb2"""                                                                                                
    cw = nmt.NmtCovarianceWorkspace()
    if not os.path.exists(cws_path):
        #raise Exception("Couldn't find an existing covariance workspace...")
        cw.compute_coupling_coefficients(fa1, fb1, fa2, fb2)
        cw.write_to(cws_path)
    else:
        cw.read_from(cws_path)
    cov = nmt.gaussian_covariance(
        cw,
        0,
        0,
        0,
        0,  # Spins of the 4 fields
        # EE, EB, BE, BB
        cla1b1,
        cla1b2,
        cla2b1,
        cla2b2,
        wa,
        wb=wb,
    ).reshape([n_ell, 1, n_ell, 1])
    return cov

# Calculate ggl-ggl gaussian covariance
def ggl_gcov_make(
    fa1, fa2, fb1, fb2, wa, wb, cla1b1, cla1b2, cla2b1, cla2b2, n_ell, cws_path
):
    """Returns the Gaussian covariance matrix fa1fa2_fb1fb2"""
    cw = nmt.NmtCovarianceWorkspace()
    if not os.path.exists(cws_path):
        raise Exception("Couldn't find an existing covariance workspace...")
        cw.compute_coupling_coefficients(fa1, fb1, fa2, fb2)
        cw.write_to(cws_path)
    else:
        cw.read_from(cws_path)
    cov = nmt.gaussian_covariance(
        cw,
        0,
        2,
        0,
        2,
        cla1b1,
        cla1b2,
        cla2b1,
        cla2b2,
        wa,
        wb=wb,
    ).reshape([n_ell, 2, n_ell, 2])
    return cov

# Calculate clustering-ggl gaussian covariance
def gclxggl_gcov_make(
    fa1, fa2, fb1, fb2, wa, wb, cla1b1, cla1b2, cla2b1, cla2b2, n_ell, cws_path
):
    """Returns the Gaussian covariance matrix fa1fa2_fb1fb2"""
    cw = nmt.NmtCovarianceWorkspace()
    if not os.path.exists(cws_path):
        raise Exception("Couldn't find an existing covariance workspace...")
        cw.compute_coupling_coefficients(fa1, fb1, fa2, fb2)
        cw.write_to(cws_path)
    else:
        cw.read_from(cws_path)
    cov = nmt.gaussian_covariance(
        cw,
        0,
        0,
        0,
        2,
        cla1b1,
        cla1b2,
        cla2b1,
        cla2b2,
        wa,
        wb=wb,
    ).reshape([n_ell, 1, n_ell, 2])
    return cov

#-----------------------------------------------
# Defining paths

data_path = "/global/cscratch1/sd/ljfaga/DESY3_data/"
theo_path = "/global/homes/l/ljfaga/3x2hs_measurements/conf/Cl_flaskv2p0_nolimber_emu_Nsource4/"
cws_path = '/global/cscratch1/sd/ljfaga/DESY3_data/cov_workspace/'
cov_path = '/global/homes/l/ljfaga/3x2hs_measurements/gcov_3x2pt/'

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

gclgcl_pairs = [(gcl_pairs[a], gcl_pairs[b])
             for a in range(len(gcl_pairs))
             for b in range(len(gcl_pairs))]

gclggl_pairs = [(gcl_pairs[a], ggl_pairs[b])
             for a in range(len(gcl_pairs))
             for b in range(len(ggl_pairs))]

gglggl_pairs = [(ggl_pairs[a], ggl_pairs[b])
             for a in range(len(ggl_pairs))
             for b in range(a, len(ggl_pairs))]

#-------------------------------------------
# Creating and calculating (or loading) workspaces
# TO DO: calculate csh and ggl workspaces if they don't exist
print('Creating and calculating (or loading) workspaces\n')

print('Shear\n')
if os.path.exists(f'{data_path}/redmagic/mcmwsp_csh_y3data_00.fits'):
    print('csh workspace already exists')
    w_csh = [
        mcm_read(f'{data_path}/redmagic/mcmwsp_csh_y3data_{i}{j}.fits')
            for i, j in csh_pairs
    ]
else:
    raise Exception("Couldn't find csh workspaces.")
#else:
#w_csh = [csh.mcm_make(field[i], field[j], bins) for i, j in pairs]

print('Clustering\n')
# Clustering
#if os.path.exists(f'{data_path}/redmagic/mcmwsp_gcl_y3data_00.fits'):
#    print('já tem gcl workspace')
#    w_gcl = nmt.NmtWorkspace()
#    w_gcl.read_from(f'{data_path}/redmagic/mcmwsp_gcl_y3data_00.fits')
#else:
w_gcl = mcm_make(conf["redmagic_mask"], gclmask, bins, odir)

print('Galaxy-galaxy lensing\n')
if os.path.exists(f'{data_path}/redmagic/mcmwsp_ggl_y3data_00.fits'):
    print('já tem ggl workspace')
    w_ggl = [mcm_read(f'{data_path}/redmagic/mcmwsp_ggl_y3data_{i}{j}.fits') 
             for i, j in ggl_pairs]
else: 
    raise Exception("Couldn't find ggl workspaces.")
#else:
#    w_ggl = [w.compute_coupling_matrix(field_i[i], cshfield[j], bins) for i ...]

#-------------------------------------------
# Loading theoretical cls
print('Loading theoretical C_ells\n')

lmax = 3 * conf["nside"]

print('Shear C_ells\n')
cl_csh = [[None for i in range(conf["nz_src"])] for j in range(conf["nz_src"])]
for ii, (i, j) in enumerate(csh_pairs):
    print(ii)
    cl = load_inputcl(i+1, j+1, lmax, idir=theo_path, cl_type='csh')
    print(cl.shape, w_csh[ii].wsp.lmax + 1, type(w_csh[ii].wsp.lmax + 1), ii, type(ii))
    cl[:, : w_csh[ii].wsp.lmax + 1] = w_csh[ii].couple_cell(cl)
    cl[:, w_csh[ii].wsp.lmax :] = cl[:, w_csh[ii].wsp.lmax][:, None]
    if i == j:
        cl += csh.pclnoise_make(cshcat[i], cshmask[i])
    cl /= np.mean(cshmask[i] * cshmask[j])
    if conf["pixwin"]:
        cl /= np.array([hp.pixwin(conf["nside"])] * 4) ** 2
    cl_csh[i][j] = cl
    if i != j:
        cl_csh[j][i] = cl
print(len(cl_csh), len(cl_csh[0]), len(cl_csh[0][0]), cl_csh[0][0][0].shape)

# Clustering
print('Clustering C_ells\n')
cl_gcl = [None for i in range(conf["nz_lns"])]
for ii, (i, j) in enumerate(gcl_pairs):
    print(ii)
    cl = load_inputcl(i+1, j+1, lmax, idir=theo_path, cl_type='gcl')
    print(cl.shape, w_gcl.wsp.lmax + 1, type(w_gcl.wsp.lmax + 1), ii, type(ii))
    cl = w_gcl.couple_cell(cl)
    # TO DO::::::: add shot noise for  clustering here
    # cl += gcl.pclnoise_make(gclcat[i], gclmask[i])
    # cl /= np.mean(gclmask[i] * gclmask[j])
    if conf["pixwin"]:
        cl /= np.array([hp.pixwin(conf["nside"])]) ** 2
    cl_gcl[ii] = cl
print(len(cl_gcl), len(cl_gcl[0]), len(cl_gcl[0][0]), cl_gcl[0][0][0].shape)

# Galaxy-galaxy lensing
print('GGL C_ells\n')
cl_ggl = [[None for i in range(conf["nz_src"])] for j in range(conf["nz_lns"])]
for ii, (i, j) in enumerate(ggl_pairs):
    cl = load_inputcl(i+1, j+1, lmax, idir=theo_path, cl_type='ggl')
    print(cl.shape, w_ggl[ii].wsp.lmax + 1, type(w_ggl[ii].wsp.lmax + 1), ii, type(ii))
    cl[:, : w_ggl[i].wsp.lmax + 1] = w_ggl[i].couple_cell(cl)
    cl[:, w_ggl[i].wsp.lmax :] = cl[:, w_ggl[i].wsp.lmax][:, None]
    if conf["pixwin"]:
        cl /= np.array([hp.pixwin(conf["nside"])] * 2) ** 2
    cl_ggl[i][j] = cl
print(len(cl_ggl), len(cl_ggl[0]))

#-------------------------------------------------
# Calculating gaussian covariance
print('Calculating gaussian covariance terms\n')

print(f'Searching for pre-existing covariance workspaces at {cws_path}')

print(f'Saving them at {cov_path} \n')

n_ell = bins.get_n_bands()

#---#---#

if block == 'gcl-gcl':

    print('Clustering-clustering block')
    
    if os.path.exists(f'{cov_path}/gcov_gcl.npz'):
        cov = np.load(f'{cov_path}/gcov_gcl.npz')
    else:
        cov = {}
    
    covpairs = gclgcl_pairs[comb][:]

    a1, a2 = covpairs[0]
    b1, b2 = covpairs[1]
    
    print(f"{a1}{a2}_{a1}{a2}")
    cov[f"{a1}{a2}_{a1}{a2}"] = gcl_gcov_make(
        gclfield[a1],
        gclfield[a2],
        gclfield[a1],
        gclfield[a2],
        w_gcl,
        w_gcl,
        cl_gcl[a1],
        cl_gcl[a1],
        cl_gcl[a2],
        cl_gcl[a2],
        n_ell,
        cws_path
    )
    print('Finished clustering-clustering')
    np.savez_compressed(f'{cov_path}/gcov_gcl.npz', **cov)
    print('Clustering-clustering block saved as gcov_gcl.npz\n\n')

#---#---#

if block == 'gcl-gcl':

    print('clustering-GGL block')
    
    if os.path.exists(f'{cov_path}/gcov_gclxggl.npz'):
        cov = np.load(f'{cov_path}/gcov_gclxggl.npz')
    else:
        cov = {}
    
    covpairs = gclggl_pairs[comb][:]

    a1, a2 = covpairs[0]
    b1, b2 = covpairs[1]

    print(f"{a1}{a2}_{b1}{b2}")
    cov[f"{a1}{a2}_{b1}{b2}"] = gclxggl_gcov_make(
        gclfield[a1],
        gclfield[a2],
        gclfield[b1],
        cshfield[b2],
        w_gcl,
        w_ggl[b],
        cl_gcl[a1], #[b1],
        cl_ggl[a2][b2],
        cl_gcl[a2], #[b1],
        cl_ggl[a1][b2],
        n_ell,
        cws_path
    )
    print('Finished clustering-GGL')
    np.savez_compressed(f'{cov_path}/gcov_gclxggl.npz', **cov)
    print(f'Clustering-GGL block saved as {cov_path}/gcov_gclxggl.npz\n\n')

#---#---#

if block == 'ggl-ggl':

    print('GGL-GGL block')
    
    if os.path.exists(f'{cov_path}/gcov_ggl.npz'):
        cov = np.load(f'{cov_path}/gcov_ggl.npz')
    else:
        cov = {}
    
    covpairs = gglggl_pairs[comb][:]

    a1, a2 = covpairs[0]
    b1, b2 = covpairs[1]

    print(f"{a1}{a2}_{b1}{b2}")
    cov[f"{a1}{a2}_{b1}{b2}"] = ggl_gcov_make(
        gclfield[a1],
        gclfield[b1],
        cshfield[a2],
        cshfield[b2],
        w_ggl[a],
        w_ggl[b],
        cl_gcl[a1],
        cl_ggl[a1][b2],
        cl_ggl[a2][b1],
        cl_csh[a2][b2],
        n_ell,
        cws_path
    )
    
    print('Finished GGL-GGL')
    np.savez_compressed(f'{cov_path}/gcov_ggl.npz', **cov)
    print(f'GGL-GGL block saved as {cov_path}/gcov_ggl.npz\n\n')

#---#---#

                                                                                                






