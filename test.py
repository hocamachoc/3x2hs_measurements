from sys import argv
from yaml import safe_load
from Cls import cls_flask
from os.path import (basename, exists)
from os import makedirs
from numpy import (loadtxt, savez_compressed)


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

if not exists(odir):
  makedirs(odir)

cls = cls_flask(iseed, ick, ell_ini, ell_end, conf['nside'], conf['nz'], 
                conf['flskdir'], conf['mskdir'], odir)
savez_compressed(f'{odir}/cls_gcl_s{iseed}_ck{ick}.npz', **cls)
