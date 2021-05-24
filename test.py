from sys import argv
from Cls import cls_flask
from os.path import basename, exists
from os import makedirs
from numpy import loadtxt, savez_compressed


id = int(argv[1])  # starts at 0
iseed, ick = id // 2 + 1, id % 2 + 1
nside, nz = 1024, 5
print(iseed, ick)

elledges = "etc/binLogN20Lmin30Lmax3000.txt"
flskdir = "/store/hcamacho/des_flasky3/maskedcats"
mskdir = "/store/hcamacho/des_flasky3/cookies"
odir = (
    "/store/hcamacho"
    + f"/cls_desFlasky3_nside{nside}_{basename(elledges).replace('.txt', '')}"
)

ell_ini, ell_end = loadtxt(elledges, unpack=True, dtype="i4")

if not exists(odir):
    makedirs(odir)

cls = cls_flask(iseed, ick, ell_ini, ell_end, nside, nz, flskdir, mskdir, odir)
savez_compressed(f"{odir}/cls_gcl_s{iseed}_ck{ick}.npz", **cls)
