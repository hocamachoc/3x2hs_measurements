"""
This module handles with post-processing of FLASK simulations
"""
import pandas as pd
import numpy as np
import healpy as hp
import fitsio as ft
import random as rd
import healpix_util as hu
import itertools as it
import multiprocessing as mp
import yaml
import os
from argparse import ArgumentParser


def cshcat_make(map_fn, nbar, sigma_e, fgoodmap_fn, fgood_thold=0.0,
                seed=None):
    """
    Generates a cosmic-shear catalog for a FLASK simulation realization.
    """
    cshcat = pd.DataFrame(columns=['ra', 'dec', 'g1', 'g2',
                                   'g1_true', 'g2_true'])

    gamma1, gamma2, ip_good, nside_map, fgoodmap = ggmap_load(
        map_fn, fgoodmap_fn, fgood_thold=fgood_thold)

    cshcat['ra'], cshcat['dec'], nc_map = ggmap_sample_positions(
        ip_good, nside_map, nbar, fgoodmap, seed=seed)
    dg1, dg2 = cshcat_samplegamma(sigma_e, nc_map.sum(), seed=seed)

    idx = [[i] * nc_map[i] for i in range(len(ip_good))]
    idx = [ip for sub in idx for ip in sub]

    cshcat['g1_true'], cshcat['g2_true'] = -gamma1[idx], gamma2[idx]
    cshcat['g1'] = -gamma1[idx] + dg1
    cshcat['g2'] = gamma2[idx] + dg2        

    return cshcat


def ggmap_load(map_fn, fgoodmap_fn, fgood_thold=0.0):
    """
    Loads a shear map from a weak-lensing FLASK simulated map (kggmap) and
    reduce it to a footprint specified by a fracgood map keeping only the shear
    components (ggmap). Note we're keeping the resolution, i.e., the resolution
    of the output ggmap is the same of the input kggmap. Note further all
    output maps are for "good", i.e., inside footprint pixels, not full-sky
    maps.
    """
    nside_map = ft.read_header(map_fn, 1)['nside']
    nside_msk = ft.read_header(fgoodmap_fn, 1)['nside']
    assert nside_map >= nside_msk 

    print("Reading", map_fn)
    fgoodmap = hp.read_map(fgoodmap_fn, verbose=False)
    ip_good, = np.where(fgoodmap > fgood_thold)
    gamma1, gamma2 = hp.read_map(map_fn, field=[1, 2], verbose=False)
    if nside_msk > nside_map:  # 'Down-grade' the ip_good accordingly
        ang_mid = hp.pix2ang(nside_map, np.arange(hp.nside2npix(nside_map)))
        mask_tmp = np.in1d(hp.ang2pix(nside_msk, *ang_mid), ip_good)
        ip_good = np.arange(hp.nside2npix(nside_map))[mask_tmp]

    return (gamma1[ip_good], gamma2[ip_good], ip_good, nside_map,
            fgoodmap[ip_good])


def ggmap_sample_positions(ip_good, nside, nbar, fgoodmap, seed=None,
                           nside_up=4):
    """
    Samples positions of galaxies inside a footprint matching an angular number
    density nbar. Currently, this do not take into account a frac-good map,
    i.e, assumes the whole area of pixels inside the footprint as covered by
    the footprint.
    """
    np.random.seed(seed)
    nbar_pix = nbar * hp.nside2pixarea(nside, degrees=True) * 3600.0
    nc_map = np.random.poisson(nbar_pix, len(ip_good))

    ra, dec = ncmap_sample_positions(nc_map, ip_good, nside, fgoodmap,
                                     nside_up=nside_up)
    
    return ra, dec, nc_map


def cshcat_samplegamma(sigma_e, ngal, seed=None):
    """
    Sample cosmic-shear dispersion
    """
    np.random.seed(seed)
    dgamma1 = sigma_e * np.random.randn(ngal)
    dgamma2 = sigma_e * np.random.randn(ngal)

    return dgamma1, dgamma2


def process_one_kggmap(iseed, ick, iz, flaskdir, outdir, nbar, sigma_e,
                       kggpref="kappa-gamma-f10"):
    """
    Processes a single FLASK seed
    """
    fgoodmap_fn = f"{flaskdir}/cookies/ck{ick+1}.fits.gz"
    map_fn = f"{flaskdir}/4096/seed{iseed+1}/{kggpref}z{iz+1}.fits"
    cat_fn = f"{outdir}/srccat_z{iz+1}_s{iseed+1}_ck{ick+1}.parquet"

    cshcat = cshcat_make(map_fn, nbar, sigma_e, fgoodmap_fn)
    cshcat.to_parquet(cat_fn, index=False)

    return cat_fn


def lnscat_load(iseed, flaskdir):
    """
    Loads a FLASK lens-catalog to a pandas data-frame
    """
    lnscat_fn = f"{flaskdir}/4096/seed{iseed+1}/lens-catalog.fits.gz"
    lnscat = ft.read(lnscat_fn, columns=['RA', 'DEC', 'galtype'])
    # lnscat = lnscat.byteswap().newbyteorder()
    lnscat = pd.DataFrame.from_records(lnscat)
    names = {'galtype': 'zbin', 'RA': 'ra', 'DEC': 'dec'}
    return lnscat.rename(columns=names)


def lnscat_addck(lnscat, flaskdir, nck=2, nz_lns=5, fgood_thold=0.0):
    """
    Adds cookie-cut information to lnscat
    """
    fgoodmap_fn = [f"{flaskdir}/cookies/ck{ick+1}_desy3_goldv2p2p1.fits.gz"
                   for ick in range(nck)]
    nside = ft.read_header(fgoodmap_fn[0], 1)['nside']
    assert nside == ft.read_header(fgoodmap_fn[1], 1)['nside']
    npix = hp.nside2npix(nside)
    
    ipix = hu.HealPix('ring', nside).eq2pix(lnscat.ra, lnscat.dec)
    ck = np.zeros(lnscat.shape[0], dtype='i4')
    for ick in range(nck):
        fgmap = hp.read_map(fgoodmap_fn[ick], verbose=False)
        ip_good = np.arange(npix)[(fgmap > fgood_thold)]
        ck[np.in1d(ipix, ip_good)] = ick + 1
    lnscat['ck'] = ck
    lnscat = lnscat[lnscat.ck > 0]
        
    return lnscat

def process_lenscat(iseed, flaskdir, outdir, nz_lns=5, nck=2):
    """
    Processes a single seed FLASK lens-catalog
    """
    lnscat = lnscat_load(iseed, flaskdir)
    lnscat = lnscat_addck(lnscat, flaskdir)

    for iz, ick in it.product(range(nz_lns), range(nck)):
        cat_fn = f"{outdir}/lnscat_z{iz+1}_s{iseed+1}_ck{ick+1}.parquet"
        lnscat[(lnscat.zbin == iz + 1) & (lnscat.ck == ick + 1)].\
            drop(columns=['zbin', 'ck']).\
            to_parquet(cat_fn, index=False)
    return


def pmap_load(map_fn, fgoodmap_fn, fgood_thold=0.0):
    """
    Loads a FLASK number counts map
    """
    nside_map = ft.read_header(map_fn, 1)['nside']
    nside_msk = ft.read_header(fgoodmap_fn, 1)['nside']
    assert nside_map >= nside_msk 

    print("Reading", map_fn)
    fgoodmap = hp.read_map(fgoodmap_fn, verbose=False)
    ip_good, = np.where(fgoodmap > fgood_thold)
    nc = hp.read_map(map_fn, verbose=False, dtype=int)
    if nside_msk > nside_map:  # 'Down-grade' the ip_good accordingly
        ang_mid = hp.pix2ang(nside_map, np.arange(hp.nside2npix(nside_map)))
        mask_tmp = np.in1d(hp.ang2pix(nside_msk, *ang_mid), ip_good)
        ip_good = np.arange(hp.nside2npix(nside_map))[mask_tmp]

    return (nc[ip_good], ip_good, nside_map, fgoodmap[ip_good])


def ncmap_sample_positions(nc_map, ip_good, nside, fgoodmap, nside_up=4):
    """
    Samples positions of galaxies inside a footprint from a number counts map
    """
    ip_good_nest = hp.ring2nest(nside, ip_good)
    ipix_nest = [rd.sample(range(ip_good_nest[i] * 4**nside_up,
                                 (ip_good_nest[i] + 1) * 4**nside_up),
                           nc_map[i])
                 for i in range(len(ip_good))]
    ipix_nest = [ip for sub in ipix_nest for ip in sub]

    hpix = hu.HealPix('nest', nside * 2**nside_up)
    ra, dec = hpix.pix2eq(ipix_nest)
    return ra, dec


def lnscat_make_frompmap(map_fn, fgoodmap_fn, fgood_thold=0.0, seed=None):
    """
    Generates a lens position catalog for a FLASK simulation realization.
    """
    lnscat = pd.DataFrame(columns=['ra', 'dec'])

    nc_map, ip_good, nside_map, fgoodmap = pmap_load(map_fn, fgoodmap_fn,
                                                     fgood_thold=fgood_thold)

    lnscat['ra'], lnscat['dec'] = ncmap_sample_positions(
        nc_map, ip_good, nside_map, fgoodmap)

    return lnscat


def process_one_pmap(iseed, ick, iz, flaskdir, outdir):
    """
    Processes a single FLASK seed p-map
    """
    fgoodmap_fn = f"{flaskdir}/cookies/ck{ick+1}.fits.gz"
    map_fn = f"{flaskdir}/4096/seed{iseed+1}/p-s{iseed+1}-f1z{iz+1}.fits"
    cat_fn = f"{outdir}/lnscat_z{iz+1}_s{iseed+1}_ck{ick+1}.parquet"

    lnscat = lnscat_make_frompmap(map_fn, fgoodmap_fn)
    lnscat.to_parquet(cat_fn, index=False)

    return cat_fn


def load_inputcl(i, j, lmax,
                 idir="/global/cscratch1/sd/faoli/flask_desy3/Cl_flaskv2p0_nolimber_emu_Nsource4"):
    """
    Loads FLASK input Cls to a numpy array
    """
    fn = f'{idir}/Y3_5x2pt_Nsource4-Cl_f10z{i+1}f10z{j+1}.dat'
    # TODO: Introduce a check here. Normally, l starts at 1 on FLASK/CLike
    # predictions, but we can introduce a check to allow for generality.
    cl = np.loadtxt(fn, usecols=1)
    cl = np.insert(cl, 0, 0)[:lmax]
    return np.array([cl] + [np.zeros_like(cl)] * 3)


def load_inputcl_y1(i, j, lmax,
                 idir="/global/cscratch1/sd/hcamacho/flask_desy1/Cl_flask"):
    """
    Loads FLASK input Cls to a numpy array
    """
    fn = f'{idir}/Y15x2pt-Cl_f2z{i+1}_f2z{j+1}'
    # TODO: Introduce a check here. Normally, l starts at 1 on FLASK/CLike
    # predictions, but we can introduce a check to allow for generality.
    cl = np.loadtxt(fn, usecols=1)
    cl = np.insert(cl, 0, 0)[:lmax]
    return np.array([cl] + [np.zeros_like(cl)] * 3)


if __name__ == "__main__":
    parser = ArgumentParser(description="FLASK post-processing")
    parser.add_argument("conf", nargs=1,
                        help="YAML configuration file")
    parser.add_argument("--iseed", type=int, required=True,
                        help="Seed index (0-based).")
    parser.add_argument("--des_release", required=True,
                        help="y1 or y3 for the moment")
    parser.add_argument("--processes", type=int,
                        default=10,
                        help="Number of processes for multiprocessing.")
    o = parser.parse_args()

    with open(o.conf[0], 'rt') as f:
        conf = yaml.safe_load(f)

    outdir = f"{conf['flaskdir']}/maskedcats"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # lnscat generation depends on if Y1 or Y3 release.
    # Also the prefix names of FLASK products depend on that.
    if o.des_release == 'y1':
        kggpref = f"kgg-s{o.iseed+1}-f2"

        args = [(o.iseed, ick, iz, conf['flaskdir'], outdir) for iz, ick
                in it.product(range(conf['nz_lns']), range(conf['nck']))]
        with mp.Pool(processes=o.processes) as pool:
            cat_fn = pool.starmap(process_one_pmap, args)

    elif o.des_release == 'y3':
        kggpref = f"kappa-gamma-f10"
        process_lenscat(o.iseed, conf['flaskdir'], outdir)

    else:
        raise NotImplementedError("Currently only y1 and y3 des_release(s) "
                                  + "are supported.")

    args = [(o.iseed, ick, iz, conf['flaskdir'], outdir, conf['neff'][iz],
             conf['sigma_e'][iz], kggpref) for iz, ick
            in it.product(range(conf['nz_src']), range(conf['nck']))]
    with mp.Pool(processes=o.processes) as pool:
        cat_fn = pool.starmap(process_one_kggmap, args)
