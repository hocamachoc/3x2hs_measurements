from fitsio import FITS
from pandas import DataFrame, concat
from healpix_util import HealPix
from healpy import npix2nside, read_map
from numpy import argwhere, isin, ones
from multiprocessing import Pool
from os.path import exists

mcalpath = 'mcalcat.fits.gz'


# def mcalcat_process(cat, zbin, mask, nside=4096):
def mcalcat_process(cat, zbin, nside, fgoodmap):
    if exists(mcalpath):
        mcat = mcalcat_read(mcalpath)
        return mcat
    cat = mcalcat_load(cat)
    # cat = mcalcat_prune(cat, fgoodmap)
    cat = mcalcat_addip(cat, nside)
    cat = mcalcat_addzbin(cat, zbin)
    cat = mcalcat_addR(cat)
    cat = mcalcat_splitzbins(cat)
    with Pool(processes=len(cat)) as pool:
        mcat = pool.map(mcal_addgg, cat)
    mcalcat_write(mcat, mcalpath)

    return mcat


def mcalcat_read(path):
    cat = []
    fits = FITS(path, 'r')
    for c in fits[1:]:
        cat.append(DataFrame.from_records(
            c.read().byteswap().newbyteorder()).set_index('coadd_objects_id'))
    fits.close()
    return cat


def mcalcat_write(cat, path):
    from os import remove
    if exists(path):
        remove(path)
    fits = FITS(path, 'rw')
    for c in cat:
        fits.write(c.reset_index().to_records(index=False))
    fits.close()


def mcal_addgg(d):
    d['g1'] = d['e1'].sub(d['e1'].mean())
    d['g2'] = d['e2'].sub(d['e2'].mean())
    d['varg'] = 0.5 * (d['e1']**2 + d['e2']**2)
    return d.drop(['e1', 'e2'], axis=1)


def mcalcat_splitzbins(d):
    zbin = sorted(d['zbin'].unique().tolist())
    return [d[d['zbin'] == z].drop('zbin', axis=1) for z in zbin]


def mcalcat_addzbin(d, zbin):
    fits = FITS(zbin)
    zbin = fits[1].read(columns=['coadd_objects_id', 'zbin_mcal'])
    fits.close()
    zbin = zbin.byteswap().newbyteorder()
    zbin = DataFrame.from_records(zbin).set_index('coadd_objects_id')
    zbin = zbin.rename(columns={"zbin_mcal": "zbin"})
    cat = concat([d, zbin], axis=1)
    return cat[cat['zbin'] > -1].dropna()


def mcalcat_addR(d):
    aux = d['R11'].add(d['R22'])
    d['R'] = aux.div(2.0)
    return d.drop(['R11', 'R22'], axis=1)


def mcalcat_prune(d, fgoodmap_fn, fgood_thold=0.0):
    fgoodmap = read_map(fgoodmap_fn, verbose=False)
    nside = npix2nside(len(fgoodmap))
    ipgood = argwhere(fgoodmap > fgood_thold).flatten()
    if f'ip{nside}' not in d:
        ip = mcalcat_addip(d, nside)[f'ip{nside}'].values
    else:
        ip = d[f'ip{nside}'].values
    m = isin(ip, ipgood)
    return d[m]


def mcalcat_addip(d, nside):
    hp = HealPix('ring', nside)
    d[f'ip{nside}'] = hp.eq2pix(d.ra.values, d.dec.values).astype('i8')
    return d.drop(['ra', 'dec'], axis=1)


def mcalcat_load(cat):
    fields = ['coadd_objects_id', 'ra', 'dec',
              'e1', 'e2', 'R11', 'R22', 'flags_select']
    fits = FITS(cat)
    cat = fits[1].read(columns=fields)
    fits.close()
    cat = cat.byteswap().newbyteorder()
    cat = DataFrame.from_records(cat).set_index('coadd_objects_id')
    # Select usable objects
    cat = cat[cat.flags_select == 0]
    # Weights
    cat['w'] = ones(cat.shape[0])
    return cat.drop(['flags_select'], axis=1)
