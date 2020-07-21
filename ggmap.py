from numpy import bincount
from healpy import nside2npix, read_map, write_map
from os.path import exists
from multiprocessing import Pool


def ggmap_process(ggmap_pref, cat, nside=4096):
	path = [f"{ggmap_pref}{ii+1}.fits" for ii in range(len(cat))]
#	with Pool(processes=len(cat)) as pool:
#		ggmap = pool.starmap(ggmap_process_single, [(path[ii], cat[ii], nside) for ii in range(len(cat))])
	ggmap = [ggmap_process_single(path[ii], cat[ii], nside) for ii in range(len(cat))]
	return ggmap


def ggmap_process_single(path, cat, nside):
	if exists(path):
		return ggmap_read_from(path)
	ggmap = ggmap_cat2ggmap(cat, nside)
	ggmap_write_to(ggmap, path)
	return ggmap


def ggmap_cat2ggmap(cat, nside):
	npix = nside2npix(nside)
	ncmap = bincount(cat[f'ip{nside}'].values, minlength=npix)
	g1 = bincount(cat[f'ip{nside}'].values, weights=cat['g1'].values, minlength=npix)
	g2 = bincount(cat[f'ip{nside}'].values, weights=cat['g2'].values, minlength=npix)
	m = (ncmap > 0)
	g1[m] = g1[m] / ncmap[m]
	g2[m] = g2[m] / ncmap[m]
	return [g1, g2]


def ggmap_write_to(gg, path):
	write_map(path, gg, overwrite=True)


def ggmap_read_from(path):
	gg = read_map(path, verbose=False, field=None)
	assert gg.shape[0] == 2
	return [gg[0], gg[1]]

