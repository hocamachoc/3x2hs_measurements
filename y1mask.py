from healpy import read_map, npix2nside, ud_grade, UNSEEN
from numpy import allclose


def y1mask_process(path, nside=4096):
    mask, fsky = y1mask_load(path, nside)
    return mask, fsky


def y1mask_load(path, nside):
    m = read_map(path, verbose=False)
    nsidein = npix2nside(len(m))
    assert nside <= nsidein
    m[(m == UNSEEN)] = 0.0
    fskyin = m.sum() / len(m)
    if nside == nsidein:
        return m, fskyin
    m = ud_grade(m, nside)
    fsky = m.sum() / len(m)
    assert allclose(fsky, fskyin, rtol=1e-5, atol=1e-8)
    return m, fsky
