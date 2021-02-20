"""Module for handling cosmic-shear maps.
"""
import pandas as pd
import numpy as np
import healpy as hp
import healpix_util as hu
import pymaster as nmt


def cat_fromflsk(cshcat_fn, nside, nonoise=False):
    """Load FLASK cosmic shear catalog and adds additional info
    """
    cshcat = pd.read_parquet(cshcat_fn)

    # Shear
    if nonoise:
        cshcat['g1'] = cshcat['g1_true']
        cshcat['g2'] = cshcat['g2_true']
    cshcat.drop(['g1_true', 'g2_true'], axis=1, inplace=True)
    
    # Angular positions
    cshcat[f'ip{nside}'] = hu.HealPix('ring', nside).eq2pix(cshcat['ra'],
                                                            cshcat['dec'])
    cshcat.drop(['ra', 'dec'], axis=1, inplace=True)

    # Extra information - here all simulated
    cshcat['R'] = np.ones(cshcat.shape[0])
    cshcat['w'] = np.ones(cshcat.shape[0])

    # Apply corrections
    cshcat['g1'] -= np.average(cshcat['g1'], weights=cshcat['w'])
    cshcat['g2'] -= np.average(cshcat['g2'], weights=cshcat['w'])

    # Shape noise
    cshcat['varg'] = 0.5 * (cshcat['g1']**2 + cshcat['g2']**2)
    
    return cshcat


def mask_make(cshcat, nside):
    """Generates a cosmic shear mask.
    """
    wgtmap = np.bincount(cshcat[f'ip{nside}'], minlength=hp.nside2npix(nside),
                         weights=cshcat['w'])
    return wgtmap


def field_make(cshcat, cshmask, purify_e=False, purify_b=False):
    """Generates NaMASTER cosmic-shear field
    """
    npix = cshmask.shape[0]
    nside = hp.npix2nside(npix)

    wg1map = np.bincount(cshcat[f'ip{nside}'], minlength=npix,
                         weights=cshcat['w'] * cshcat['g1'])
    wg2map = np.bincount(cshcat[f'ip{nside}'], minlength=npix,
                         weights=cshcat['w'] * cshcat['g2'])

    Rbias_mean = ((cshcat['R'] * cshcat['w']).sum() / cshcat['w'].sum())

    ip_good = (cshmask > 0.0)
    q = np.zeros_like(wg1map)
    q[ip_good] = wg1map[ip_good] / (cshmask[ip_good] * Rbias_mean)
    u = np.zeros_like(wg2map)
    u[ip_good] = wg2map[ip_good] / (cshmask[ip_good] * Rbias_mean)

    return nmt.NmtField(cshmask, [-q, u])


def mcm_make(cshfld_0, cshfld_1, bins):
    """ Return MCM workspace for shear - shear
    """
    w = nmt.NmtWorkspace()
    w.compute_coupling_matrix(cshfld_0, cshfld_1, bins)
    return w


def pclnoise_make(cshcat, cshmask):
    """Generates NaMASTER PCLs for shape noise
    """
    npix = cshmask.shape[0]
    nside = hp.npix2nside(npix)

    w2varg_pixmean = (cshcat['w']**2 * cshcat['varg']).sum() / npix
    Rbias_mean = ((cshcat['R'] * cshcat['w']).sum() / cshcat['w'].sum())

    pix_area = 4.0 * np.pi / npix
    pcl_noise = np.full(3 * nside, pix_area * w2varg_pixmean / Rbias_mean**2)
    return np.array([pcl_noise, np.zeros(pcl_noise.shape[0]),
                     np.zeros(pcl_noise.shape[0]), pcl_noise])
