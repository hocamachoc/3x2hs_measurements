from os.path import (exists, getsize)
from healpix import (read_map, write_map, npix2nside, ud_grade)


def make_y3mask(ick, nside, out_dir, nside_in=4096):
    """ Returns the Y3 mask
    """
    assert nside <= nside_in
    path = f"{out_dir}/mask_ck{ick}_nside{nside}.fits"
    if exists(path) and getsize(path) > 0:
        msk = read_map(path, verbose=False)
        return msk
    msk = process_y3mask(ick, nside, nside_in)
    write_map(path, msk)
    return msk


def process_y3mask(y3mask_path, ick, nside, nside_in):
    """ Get the original Y3 mask DES product and process it to our desired format
    """
    msk = read_map(y3mask_path)
    assert nside == npix2nside(len(msk))
    if nside == nside_in:
        return msk
    msk = ud_grade(msk, nside)
    return msk
