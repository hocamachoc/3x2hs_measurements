from os.path import (exists, getsize)
from healpy import (read_map, write_map, npix2nside, ud_grade)


def make_y3mask(ick, nside, out_dir, maskdir):
    """ Returns the Y3 mask
    """
    path = f"{out_dir}/mask_gcl_ck{ick}.fits"
    if exists(path) and getsize(path) > 0:
        msk = read_map(path, verbose=False)
        return msk
    y3mask_path = f"{maskdir}/ck{ick}_desy3_goldv2p2p1.fits.gz"
    msk = process_y3mask(y3mask_path, ick, nside)
    write_map(path, msk)
    return msk


def process_y3mask(y3mask_path, ick, nside):
    """ Get the original Y3 mask DES product and process it to our desired format
    """
    msk = read_map(y3mask_path)
    nside_in = npix2nside(len(msk))
    assert nside <= nside_in
    if nside == nside_in:
        return msk
    msk = ud_grade(msk, nside)
    return msk
