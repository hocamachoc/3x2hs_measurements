"""
Module for handling galaxy clustering fields
"""
import pymaster as nmt
import healpix_util as hu
import healpy as hp
import numpy as np
import pandas as pd
import os


def cat_fromflsk(gclcat_fn, nside):
    """
    Returns a galaxy-clustering catalog from FLASK
    """
    print("Reading", gclcat_fn)
    gclcat = pd.read_parquet(gclcat_fn)
    gclcat = _add_ipix(gclcat, nside)
    # Extra information - here all simulated
    gclcat["w"] = np.ones(gclcat.shape[0])
    return gclcat


def _add_ipix(gclcat, nside):
    """
    Adds IP'NSIDE' to a dataframe originally containing {ra, dec}
    """
    hp = hu.HealPix("ring", nside)
    gclcat[f"ip{nside}"] = hp.eq2pix(gclcat.ra.values, gclcat.dec.values)
    return gclcat.drop(["ra", "dec"], axis=1)


def mask_make(mask_fn, ick, nside, cachedir):
    """
    Galaxy clustering mask
    """
    path = f"{cachedir}/mask_gcl_ck{ick}.fits"
    if os.path.exists(path) and os.path.getsize(path) > 0:
        msk = hp.read_map(path, verbose=False)
        return msk
    msk = _mask_process(mask_fn, nside)
    hp.write_map(path, msk)
    return msk


def _mask_process(mask_fn, nside):
    """
    Processes original mass
    """
    msk = hp.read_map(mask_fn)
    nside_in = hp.npix2nside(msk.shape[0])
    assert nside <= nside_in
    if nside == nside_in:
        return msk
    msk = hp.ud_grade(msk, nside)
    return msk


def mcm_make(mask, ick, b, cachedir):
    """
    Returns MCM workspace for position-position (clustering)
    """
    path = f"{cachedir}/mcmwsp_gcl_ck{ick}.fits"
    w = nmt.NmtWorkspace()
    if os.path.exists(path) and os.path.getsize(path) > 0:
        w.read_from(path)
        return w
    f1 = nmt.NmtField(mask, [mask])
    w.compute_coupling_matrix(f1, f1, b)
    w.write_to(path)
    return w


def field_make(gclcat, gclmask, nside, save_maps=False, maps_prefix=""):
    """
    Generates NaMASTER galaxy-clustering field
    Also returns the total number of objects
    """
    ipgood = gclmask > 0
    ncmap = np.bincount(
        gclcat[f"ip{nside}"], weights=gclcat["w"], minlength=len(gclmask)
    )
    dmap = np.full(gclmask.shape[0], 0.0)
    dmap[ipgood] = (
        ncmap[ipgood]
        / gclmask[ipgood]
        * gclmask[ipgood].sum()
        / ncmap[ipgood].sum()
        - 1
    )

    if save_maps:
        hp.write_map(f"{maps_prefix}_ncmap.fits", ncmap, overwrite=True)
        hp.write_map(f"{maps_prefix}_dmap.fits", dmap, overwrite=True)

    return nmt.NmtField(gclmask, [dmap]), ncmap[ipgood].sum()


def pclnoise_make(fsky, nobj, nside):
    """
    Generates NaMASTER PCLs for shot noise
    """
    ndens = fsky * 4.0 * np.pi / nobj
    nl = fsky * ndens
    return np.array([np.full(3 * nside, nl)])
