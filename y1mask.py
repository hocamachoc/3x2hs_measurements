from healpy import read_map, npix2nside, ud_grade, UNSEEN


def y1mask_process(path, nside=4096):
	mask, fsky = y1mask_load(path, nside)
	y1mask_fig(mask, fsky, '/home/hcamacho/public_html/test.png')	
	return mask, fsky


def y1mask_load(path, nside):
	m = read_map(path, verbose=False)
	nsidein = npix2nside(len(m))
	assert nside <= nsidein
	m[(m==UNSEEN)] = 0.0
	fskyin = m.sum() / len(m)
	if nside == nsidein:
		return m, fskyin
	m = ud_grade(m, nside)
	fsky = m.sum() / len(m)
	assert allclose(fsky, fskyin, rtol=1e-5, atol=1e-8)
	return m, fsky


def y1mask_fig(mask, fsky, ofn):
	from healpy import mollview
	from matplotlib.pyplot import savefig
	mollview(mask, title=r'$f_{\rm sky} =' + f'{fsky*100:.3f}' + r'\times 10^{-2}$')
	savefig(ofn)

