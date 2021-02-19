"""Module handling post-processing of FLASK simulations
"""
import pandas as pd
import numpy as np
import healpy as hp
import fitsio as ft
import time as tm
import random as rd
import healpix_util as hu
import itertools as it
import multiprocessing as mp
import os
from argparse import ArgumentParser
import sys


def cshcat_make(map_fn, nbar, sigma_e, fgoodmap_fn, fgood_thold=0.0,
                seed=None):
    """Generates a cosmic-shear catalog for FLASK simulation realization.
    """
    cshcat = pd.DataFrame(columns=['RA', 'DEC', 'GAMMA1', 'GAMMA2',
                                   'GAMMA1_TRUTH', 'GAMMA2_TRUTH'])

    gamma1, gamma2, ip_good, nside_map, fgoodmap = ggmap_load(
        map_fn, fgoodmap_fn, fgood_thold=fgood_thold)
    cshcat['RA'], cshcat['DEC'], nc_map = ggmap_sample_positions(
        ip_good, nside_map, nbar, fgoodmap, seed=seed)
    (cshcat['GAMMA1_TRUTH'], cshcat['GAMMA2_TRUTH'],
     cshcat['GAMMA1'], cshcat['GAMMA2']) = cshcat_samplegamma(
        sigma_e, gamma1, gamma2, nc_map, ip_good, nside_map, seed=seed)

    return cshcat


def ggmap_load(map_fn, fgoodmap_fn, fgood_thold=0.0):
    """Loads a shear map from a weak-lensing FLASK map (kggmap) and reduce it
    to a footprint specified by a fracgood map.
    Note the resolution of output ggmap is the same of the input kggmap.
    """
    print("> Load data")
    sec = tm.time()
    nside_map = ft.read_header(map_fn, 1)['nside']
    nside_msk = ft.read_header(fgoodmap_fn, 1)['nside']
    fgoodmap = hp.read_map(fgoodmap_fn, verbose=False)
    ip_good, = np.where(fgoodmap > fgood_thold)
    gamma1, gamma2 = hp.read_map(map_fn, field=[1, 2], verbose=False)
    assert nside_map >= nside_msk
    if nside_msk != nside_map:
        ang_mid = hp.pix2ang(nside_map, np.arange(hp.nside2npix(nside_map)))
        mask_tmp = np.in1d(hp.ang2pix(nside_msk, *ang_mid), ip_good)
        ip_good = np.arange(hp.nside2npix(nside_map))[mask_tmp]
    print(tm.time() - sec)
    return (gamma1[ip_good], gamma2[ip_good], ip_good, nside_map,
            fgoodmap[ip_good])


def ggmap_sample_positions(ip_good, nside, nbar, fgoodmap, seed=None,
                           nside_up=4):
    """Samples the positions of galaxies
    """
    np.random.seed(seed)
    nbar_pix = nbar * hp.nside2pixarea(nside, degrees=True) * 3600.0
    nc_map = np.random.poisson(nbar_pix, len(ip_good))
    
    print("Sample Positions")
    print("nbar (_pix) =", nbar, nbar_pix, nc_map.mean())
    print("Nsrc =", nc_map.sum())
    print("Aeff =", fgoodmap.shape[0] *
          hp.nside2pixarea(nside, degrees=True))
    print("  neff =", nc_map.sum() / fgoodmap.shape[0] /
          hp.nside2pixarea(nside, degrees=True) / 3600.0)
    print("  Aeff =", fgoodmap.sum() * hp.nside2pixarea(nside, degrees=True))
    print("  neff =", nc_map.sum() / fgoodmap.sum() /
          hp.nside2pixarea(nside, degrees=True) / 3600.0)

    ip = [[ip_good[i]] * nc_map[i] for i in range(len(ip_good))]
    ip = [i for sub in ip for i in sub]
    return *hu.HealPix('ring', nside).pix2eq(ip), nc_map

    # s = tm.time()
    # ip_good_nest = hp.ring2nest(nside, ip_good)
    # print(tm.time() - s)
    # s = tm.time()
    # ipix_nest = [rd.sample(range(ip_good_nest[i] * 4**nside_up,
    #                              (ip_good_nest[i] + 1) * 4**nside_up),
    #                        nc_map[i])
    #              for i in range(len(ip_good))]
    # ipix_nest = [ip for sub in ipix_nest for ip in sub]
    # print(tm.time() - s)

    # hpix = hu.HealPix('nest', nside * 2**nside_up)
    # return *hpix.pix2eq(ipix_nest), nc_map


def ggmap_sample_positions_mp(ip_good, nside, nbar, fgoodmap, seed=None,
                              nside_up=4):
    """Samples the positions of galaxies
    """
    np.random.seed(seed)
    nbar_pix = nbar * hp.nside2pixarea(nside, degrees=True) * 3600.0
    nc_map = np.random.poisson(nbar_pix, len(ip_good))
    
    print("> Sample Positions")
    print("  nbar (_pix) =", nbar, nbar_pix, nc_map.mean())
    print("  Nsrc =", nc_map.sum())
    print("  Aeff =", fgoodmap.shape[0] *
          hp.nside2pixarea(nside, degrees=True))
    print("  neff =", nc_map.sum() / fgoodmap.shape[0] /
          hp.nside2pixarea(nside, degrees=True) / 3600.0)
    print("  Aeff =", fgoodmap.sum() * hp.nside2pixarea(nside, degrees=True))
    print("  neff =", nc_map.sum() / fgoodmap.sum() /
          hp.nside2pixarea(nside, degrees=True) / 3600.0)

    s = tm.time()
    ip_good_nest = hp.ring2nest(nside, ip_good)
    print(tm.time() - s)
    s = tm.time()
    # ipix_nest = [rd.sample(range(ip_good_nest[i] * 4**nside_up,
    #                              (ip_good_nest[i] + 1) * 4**nside_up),
    #                        nc_map[i])
    #              for i in range(len(ip_good))]
    with mp.Pool() as pool:
        ipix_nest = pool.starmap(rd.sample,
                                 [(range(ip_good_nest[i] * 4**nside_up,
                                         (ip_good_nest[i] + 1) * 4**nside_up),
                                   nc_map[i])
                                  for i in range(len(ip_good))])
    ipix_nest = [ip for sub in ipix_nest for ip in sub]
    print(tm.time() - s)

    hpix = hu.HealPix('nest', nside * 2**nside_up)
    return *hpix.pix2eq(ipix_nest), nc_map


def cshcat_samplegamma(sigma_e, gamma1, gamma2, nc_map, ip_good, nside,
                       seed=None):
    """Sample shear dispersion
    """
    # Constructs a list with pixel idxs at map resolution for the
    # source galaxies. Note the idxs on gamma{1,2} and ip_good are not over all
    # map indexes, but on the masked region.
    # ipix = [[ip_good[i]] * nc_map[i] for i in range(len(ip_good))]
    ipix = [[i] * nc_map[i] for i in range(len(ip_good))]
    ipix = [ip for sub in ipix for ip in sub]
    ngal = nc_map.sum()

    np.random.seed(seed)
    dgamma1 = sigma_e * np.random.randn(ngal)
    dgamma2 = sigma_e * np.random.randn(ngal)

    return (gamma1[ipix], gamma2[ipix],
            gamma1[ipix] + dgamma1, gamma2[ipix] + dgamma2)


def process_one_kggmap(iseed, ick, iz, flaskdir, outdir, nbar, sigma_e):
    """Process a single seed
    """
    fgoodmap_fn = f"{flaskdir}/cookies/ck{ick+1}_desy3_goldv2p2p1.fits.gz"
    map_fn = f"{flaskdir}/4096/seed{iseed+1}/kappa-gamma-f10z{iz+1}.fits"
    cat_fn = f"{outdir}/srccat_z{iz+1}_s{iseed+1}_ck{ick+1}.parquet"

    cshcat = cshcat_make(map_fn, nbar, sigma_e, fgoodmap_fn)
    cshcat.to_parquet(cat_fn, index=False)

    return cat_fn
        
if __name__ == "__main__":
    sys.stdout.flush()
    
    parser = ArgumentParser()
    parser.add_argument("--iseed", type=int, required=True,
                        help="Seed index (0-based).")
    parser.add_argument("--flaskdir",
                        default="/global/cscratch1/sd/faoli/flask_desy3",
                        help="FLASK root directory.")
    parser.add_argument("--outdir",
                        default="/global/cscratch1/sd/hcamacho/flask_desy3/maskedcats",
                        help="Output directory.")
    o = parser.parse_args()

    # Y3, See cdcvs.fnal.gov/redmine/projects/des-theory/wiki/Y3_Analysis_Choices
    # May need some cross-check from Troxel.
    nbar = [1.476, 1.479, 1.484, 1.461]
    sigma_e = [0.247, 0.266, 0.263, 0.314]

    if not os.path.exists(o.outdir):
        os.makedirs(o.outdir)
    
    args = [(o.iseed, ick, iz, o.flaskdir, o.outdir, nbar[iz], sigma_e[iz])
            for iz, ick in it.product(range(len(nbar)), range(2))]
    print(args)
    with mp.Pool(processes=len(args)) as pool:
        cat_fn = pool.starmap(process_one_kggmap, args)
            
