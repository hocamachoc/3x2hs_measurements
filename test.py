from sys import argv
from yaml import safe_load
from Cls import cls_flask
from os.path import basename, exists
from os import makedirs
from numpy import loadtxt, savez_compressed


<<<<<<< HEAD
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
=======
conf = argv[1]
real_id = int(argv[2])   # Realization ID. Starts at 0 

iseed, ick = real_id // 2 + 1, real_id % 2 + 1
print(iseed, ick)

with open(conf, 'rt') as f:
    conf = safe_load(f)
odir = f"{conf['odir']}/nside{conf['nside']}"\
       f"_{basename(conf['elledges']).replace('.txt', '')}"
print(conf, odir)

ell_ini, ell_end = loadtxt(conf['elledges'], unpack=True, dtype='i4')
# TODO: Check ell_ini < 3 * nside
>>>>>>> y1ggl

if not exists(odir):
    makedirs(odir)

<<<<<<< HEAD
cls = cls_flask(iseed, ick, ell_ini, ell_end, nside, nz, flskdir, mskdir, odir)
print("Writing", f"{odir}/cls_gcl_s{iseed}_ck{ick}.npz")
savez_compressed(f"{odir}/cls_gcl_s{iseed}_ck{ick}.npz", **cls)
=======
cls = cls_flask(iseed, ick, ell_ini, ell_end, conf['nside'], conf['nz'], 
                conf['flskdir'], conf['mskdir'], odir)
savez_compressed(f'{odir}/cls_gcl_s{iseed}_ck{ick}.npz', **cls)
>>>>>>> y1ggl
