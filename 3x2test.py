import sys
import os
import yaml
import numpy as np
import healpy as hp
import pymaster as nmt
import csh
import gcl
from functools import partial
import itertools as it

print = partial(print, flush=True)  # For the impatient people :).

# Configuration file
conf = sys.argv[1]
with open(conf, "rt") as f:
    conf = yaml.safe_load(f)
odir = (
    f"{conf['odir']}/nside{conf['nside']}"
    f"_{os.path.basename(conf['elledges']).replace('.txt', '')}"
)

# Create dirs (if needed)
if not os.path.exists(odir):
    os.makedirs(odir)

if conf["type"] == "y3data":
    # Prepare cosmic-shear stuff
    cshcat_full = csh.cat_fromy3data(conf["metacal"], conf["nside"])
    cshcat = [
        cshcat_full.loc[cshcat_full["bin_number"] == iz]
        for iz in range(conf["nz_src"])
    ]

    # Prepare galaxy-clustering stuff
    gclmask = gcl.mask_make_y3data(conf["redmagic_mask"], conf["nside"])
    fsky = gclmask.mean()
    gclfield, nobj = [], []
    gclcat_full = gcl.cat_fromy3data(conf["redmagic"], conf["nside"])
    for iz in range(conf["nz_lns"]):
        gclcat = gclcat_full.loc[gclcat_full["bin_number"] == iz + 1]

        f, n = gcl.field_make(
            gclcat,
            gclmask,
            conf["nside"],
            save_maps=conf["save_maps"],
            maps_prefix=f"{odir}/zbin{iz}",
        )
        gclfield.append(f)
        nobj.append(n)

    # Output filename(s)
    ofn_csh = f"{odir}/cls_csh_y3data.npz"
    ofn_ggl = f"{odir}/cls_ggl_y3data.npz"
    ofn_gcl = f"{odir}/cls_gcl_y3data.npz"
elif conf["type"] == "flask":
    # We might want to test with true tabels
    if conf["nonoise"]:
        odir += "_nonoise"

    iseed, ick = int(sys.argv[2]), int(sys.argv[3])
    print("GGL", iseed, ick)

    # Prepare cosmic-shear stuff
    cshcat = [
        f"{conf['flaskdir']}/maskedcats"
        + f"/srccat_z{iz+1}_s{iseed}_ck{ick}.parquet"
        for iz in range(conf["nz_src"])
    ]
    cshcat = [
        csh.cat_fromflsk(fn, conf["nside"], conf["nonoise"]) for fn in cshcat
    ]

    # Prepare galaxy-clustering stuff
    gclmask = f"{conf['flaskdir']}/cookies/ck{ick}.fits.gz"
    gclmask = gcl.mask_make_flask(gclmask, ick, conf["nside"], odir)
    fsky = gclmask.mean()
    gclfield, nobj = [], []
    for iz in range(conf["nz_lns"]):
        gclcat = (
            f"{conf['flaskdir']}/maskedcats"
            + f"/lnscat_z{iz+1}_s{iseed}_ck{ick}.parquet"
        )
        gclcat = gcl.cat_fromflsk(gclcat, conf["nside"])

        f, n = gcl.field_make(
            gclcat,
            gclmask,
            conf["nside"],
            save_maps=conf["save_maps"],
            maps_prefix=f"{odir}/zbin{iz}",
        )
        gclfield.append(f)
        nobj.append(n)

    # Output filename(s)
    ofn_csh = f"{odir}/cls_csh_s{iseed}_ck{ick}.npz"
    ofn_ggl = f"{odir}/cls_ggl_s{iseed}_ck{ick}.npz"
    ofn_gcl = f"{odir}/cls_gcl_s{iseed}_ck{ick}.npz"
else:
    raise NotImplementedError(
        f"Computation type {conf['type']} not" + " implemented"
    )

# Bandpower binning - always from 0 to 3 * nside.
elledges = np.loadtxt(conf["elledges"], dtype=int)
elledges = elledges[(elledges <= 3 * conf["nside"])]
if elledges[0] > 0:
    elledges = np.insert(elledges, 0, 0)
if elledges[-1] < 3 * conf["nside"]:
    elledges = np.append(elledges, 3 * conf["nside"])
bins = nmt.NmtBin.from_edges(elledges[:-1], elledges[1:])

# CSH and GGL MCM are computed for each tomo bin pair
w = nmt.NmtWorkspace()
# GCL MCM is common to all tomo bin pair
if conf["type"] == "y3data":
    w_gcl = gcl.mcm_make(conf["redmagic_mask"], gclmask, bins, odir)
elif conf["type"] == "flask":
    w_gcl = gcl.mcm_make(
        f"{conf['flaskdir']}/cookies/ck{ick}.fits.gz", gclmask, bins, odir
    )

cls_csh = {"ell_eff": bins.get_effective_ells()}
cls_ggl = {"ell_eff": bins.get_effective_ells()}
cls_gcl = {"ell_eff": bins.get_effective_ells()}
cls_gcl["bpwrwin"] = w_gcl.get_bandpower_windows()

# CSH fields
cshfield = []
for j in range(conf["nz_src"]):
    cshcat_j = cshcat[j]
    cshmask_j = csh.mask_make(cshcat_j, conf["nside"])
    cshfield.append(
        csh.field_make(
            cshcat_j,
            cshmask_j,
            save_maps=conf["save_maps"],
            maps_prefix=f"{odir}/zbin{j}",
        )
    )

# CSH
for i, j in it.combinations_with_replacement(range(conf["nz_src"]), 2):
    w = csh.mcm_make(cshfield[i], cshfield[j], bins)
    cls_csh[f"bpwrwin_{i}{j}"] = w.get_bandpower_windows()
    cls_coup = nmt.compute_coupled_cell(cshfield[i], cshfield[j])
    if conf["pixwin"]:
        cls_coup /= np.array([hp.pixwin(conf["nside"])] * 4) ** 2
    cls_csh[f"pcl_{i}{j}"] = cls_coup
    cls_csh[f"cl_{i}{j}"] = w.decouple_cell(cls_coup)
    if i == j:
        nls_coup = csh.pclnoise_make(cshcat[i], cshfield[i].get_mask())
        if conf["pixwin"]:
            nls_coup /= np.array([hp.pixwin(conf["nside"])] * 4) ** 2
        cls_csh[f"pnl_{i}"] = nls_coup
        cls_csh[f"nl_{i}"] = w.decouple_cell(nls_coup)
    # Saving the MCM
    if conf["type"] == "flask" and iseed != 0 and ick != 1:
        continue
    w.write_to(
        ofn_csh.replace("cls_", "mcmwsp_").replace(".npz", f"_{i}{j}.fits")
    )

field_i = gclfield
for i in range(conf["nz_lns"]):
    # GCL auto-correlations
    cls_gcl[f"pcl_{i}{i}"] = nmt.compute_coupled_cell(gclfield[i], gclfield[i])
    cls_gcl[f"pnl_{i}"] = gcl.pclnoise_make(fsky, nobj[i], conf["nside"])
    if conf["pixwin"]:
        cls_gcl[f"pcl_{i}{i}"] /= np.array([hp.pixwin(conf["nside"])]) ** 2
        cls_gcl[f"pnl_{i}"] /= np.array([hp.pixwin(conf["nside"])]) ** 2
    cls_gcl[f"cl_{i}{i}"] = w_gcl.decouple_cell(cls_gcl[f"pcl_{i}{i}"])
    cls_gcl[f"nl_{i}"] = w_gcl.decouple_cell(cls_gcl[f"pnl_{i}"])
    for j in range(conf["nz_src"]):
        # GGL
        w.compute_coupling_matrix(field_i[i], cshfield[j], bins)
        cls_ggl[f"bpwrwin_{i}{j}"] = w.get_bandpower_windows()
        cls_coup = nmt.compute_coupled_cell(field_i[i], cshfield[j])
        if conf["pixwin"]:
            cls_coup /= np.array([hp.pixwin(conf["nside"])] * 2) ** 2
        cls_ggl[f"pcl_{i}{j}"] = cls_coup
        cls_ggl[f"cl_{i}{j}"] = w.decouple_cell(cls_coup)
        # Saving the MCM
        if conf["type"] == "flask" and iseed != 0 and ick != 1:
            continue
        w.write_to(
            ofn_ggl.replace("cls_", "mcmwsp_").replace(".npz", f"_{i}{j}.fits")
        )

# GCL cross-correlations
if conf["compute_cross"]:
    for i, j in it.combinations(range(conf["nz_lns"]), 2):
        cls_gcl[f"pcl_{i}{j}"] = nmt.compute_coupled_cell(
            gclfield[i], gclfield[j]
        )
        if conf["pixwin"]:
            cls_gcl[f"pcl_{i}{j}"] /= np.array([hp.pixwin(conf["nside"])]) ** 2
        cls_gcl[f"cl_{i}{j}"] = w.decouple_cell(cls_gcl[f"pcl_{i}{j}"])

print("Writing", ofn_csh)
np.savez_compressed(ofn_csh, **cls_csh)
print("Writing", ofn_ggl)
np.savez_compressed(ofn_ggl, **cls_ggl)
print("Writing", ofn_gcl)
np.savez_compressed(ofn_gcl, **cls_gcl)
