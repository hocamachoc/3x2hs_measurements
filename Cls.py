import y3flask
import y3mask
import mcm
from numpy import bincount, pi, full
from pymaster import NmtBin, compute_coupled_cell, NmtField


def cls_flask(
    iseed, ick, ell_ini, ell_end, nside, nz, flaskdir, maskdir, out_dir
):
    """Returns Cls for a FLASK realization {ISEED, ICK}"""
    b = NmtBin.from_edges(ell_ini, ell_end)
    cat = y3flask.make_poscat(iseed, ick, nside, nz, flaskdir)
    msk = y3mask.make_y3mask(ick, nside, out_dir, maskdir)
    ipgood = msk > 0
    w = mcm.make_pos(msk, ick, b, out_dir)
    cls = {"ell": b.get_effective_ells()}
    for zi in range(5):
        (
            cls[f"pcl_{zi}"],
            cls[f"pnl_{zi}"],
            cls[f"cl_{zi}"],
            cls[f"nl_{zi}"],
        ) = _get_cl(zi, cat, msk, w, nside, ipgood)
    return cls


def _get_cl(iz, cat, msk, w, nside, ipgood, wgt=False):
    """Returns auto-correlation Cls for the IZ ZBIN
    output: PCL, PNL, CL, NL
    """
    if wgt:
        nc = bincount(
            cat.loc[cat.ZBIN == iz, f"IP{nside}"],
            weights=cat.loc[cat.ZBIN == iz, "WEIGHT"],
            minlength=len(msk),
        )
    else:
        nc = bincount(
            cat.loc[cat.ZBIN == iz, f"IP{nside}"], minlength=len(msk)
        )
    f = _gc_field(nc, msk, ipgood)
    cl = compute_coupled_cell(f, f)
    ndens = msk.mean() * 4.0 * pi / nc[ipgood].sum()
    print(iz, ndens)
    nl = msk.mean() * ndens
    nl = [full(w.wsp.lmax + 1, nl)]
    nl[0][0] = 0.0
    return cl[0], nl[0], w.decouple_cell(cl)[0], w.decouple_cell(nl)[0]


def _gc_field(nc, msk, ipgood):
    """Returns over-density NmtField from number counts map"""
    dmap = full(len(msk), 0.0)
    dmap[ipgood] = (
        nc[ipgood] / msk[ipgood] * msk[ipgood].sum() / nc[ipgood].sum() - 1
    )
    return NmtField(msk, [dmap])
