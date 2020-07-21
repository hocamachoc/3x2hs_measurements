from numpy import bincount
from healpy import nside2npix
from multiprocessing import Pool


def ggmap_process(cat, nside=4096):
#	with Pool(processes=len(cat)) as pool:
#		ggmap = pool.starmap(ggmap_cat2ggmap, [(c, nside) for c in cat])
	ggmap = [ggmap_cat2ggmap(c, nside) for c in cat]
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

