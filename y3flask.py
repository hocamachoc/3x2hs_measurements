from pandas import (read_parquet, concat)
from healpix_util import HealPix


def make_poscat(iseed, ick, nside, nz, flaskdir):
    """ Returns a FLASK realization position catalog {PIXEL, ZBIN}
    """
    poscat = join_posmaskedcats(iseed, ick, nz, flaskdir)
    poscat = add_ipix_poscat(poscat, nside)
    return poscat


def join_posmaskedcats(iseed, ick, nz, flaskdir):
    """ Returns joined posmaskedcats originaly separated by ZBIN
    """
    fn = [f"{flaskdir}/lnscat_z{iz+1}_s{iseed}_ck{ick}.parquet"
          for iz in range(nz)]
    poscat = [read_parquet(f) for f in fn]
    for zi in range(nz):
        poscat[zi]['ZBIN'] = zi
    return concat(poscat)


def add_ipix_poscat(poscat, nside):
    """ Add IP'NSIDE' to a dataframe originally containing {RA, DEC}
    """
    hp = HealPix('ring', nside)
    poscat[f"IP{nside}"] = hp.eq2pix(poscat.ra.values, poscat.dec.values)
    return poscat.drop(['ra', 'dec'], axis=1)
