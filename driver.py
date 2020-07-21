from y1mask import y1mask_process
from mcalcat import mcalcat_process
from ggmap import ggmap_process
from mcm import mcm_process
from pymaster import NmtBin, NmtField, compute_coupled_cell
from numpy import loadtxt
from os.path import basename
from itertools import combinations_with_replacement, product
from pandas import DataFrame, MultiIndex
from sys import argv, exit


def main(argv):
	nside = 4096
	purify_e = False
	purify_b = False
	mask_path = "datay1/DES_Y1A1_3x2pt_redMaGiC_MASK_HPIX4096RING.fits"
	cat_path = "datay1/mcal-y1a1-combined-riz-unblind-v4-matched.fits"
	zbin_path = "datay1/y1_source_redshift_binning_v1.fits"
	bin_path = "datay1/bin_log_N20_lmin30_lmax3000.dat"

	print(f"Processing mask: {mask_path}")
	mask, fsky = y1mask_process(mask_path) 
	print(f"fsky = {fsky}, 1/fsky = {1/fsky}")

	print(f"Processing catalog: {cat_path}")
	cat = mcalcat_process(cat_path, zbin_path, mask)

	ggmap_pref = f"datay1/{basename(cat_path).split('.')[0]}_nside{nside}_ggmap_zbin"
	print(f"Processing ggmaps: {ggmap_pref}")
	ggmap = ggmap_process(ggmap_pref, cat)

	print(f"Processing binning: {bin_path}")	
	bpws, ells, weights = loadtxt(bin_path, unpack=True)
	b = NmtBin(nside=nside, bpws=bpws, ells=ells, weights=weights)

	mcm_path = f"datay1/mcm_{basename(mask_path).split('.')[0]}_{basename(bin_path).split('.')[0]}.fits"
	print(f"Processing MCM {mcm_path}")	
	w = mcm_process(mcm_path, mask, b, purify_e, purify_b)

	print(f"Processing fields")
	cols = list(product([f'{zi+1}{zj+1}' for zi, zj in combinations_with_replacement(range(len(ggmap)), 2)], ['EE', 'EB', 'BE', 'BB']))
	cols = MultiIndex.from_tuples(cols)
	f2 = [NmtField(mask, gg, purify_e=purify_e, purify_b=purify_b) for gg in ggmap]
	cl = [compute_coupled_cell(f2[zi], f2[zj]) for zi, zj in combinations_with_replacement(range(len(ggmap)), 2)]
	df = DataFrame(columns=cols)
	i = 0
	for zi, zj in combinations_with_replacement(range(len(ggmap)), 2):
		df[(f'{zi}{zj}', 'EE')] = cl[i][0]
		df[(f'{zi}{zj}', 'EB')] = cl[i][1]
		df[(f'{zi}{zj}', 'BE')] = cl[i][2]
		df[(f'{zi}{zj}', 'BB')] = cl[i][3]
		i += 1
	df.to_csv(f"datay1/coupledcls_{basename(mask_path).split('.')[0]}_{basename(bin_path).split('.')[0]}.csv.gz")
	cl = [w.decouple_cell(c) for c in cl]
	df = DataFrame(index=b.get_effective_ells(), columns=cols)
	i = 0
	for zi, zj in combinations_with_replacement(range(len(ggmap)), 2):
		df[(f'{zi}{zj}', 'EE')] = cl[i][0]
		df[(f'{zi}{zj}', 'EB')] = cl[i][1]
		df[(f'{zi}{zj}', 'BE')] = cl[i][2]
		df[(f'{zi}{zj}', 'BB')] = cl[i][3]
		i += 1
	df.to_csv(f"datay1/cls_{basename(mask_path).split('.')[0]}_{basename(bin_path).split('.')[0]}.csv.gz")

	print(f"Processing noise TODO")
	
	return 0


if __name__ == "__main__":
	exit(main(argv))

