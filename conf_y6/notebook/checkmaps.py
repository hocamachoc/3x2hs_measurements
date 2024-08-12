import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits

theo_path = '../clsdesy6/'
map_path = '/scratch/eubd/lucas.faga/tmp/tmp.1yX4PUTgWP/4096/seed1'

#m = fits.open(f'{map_path}/map-f1z1.fits')
m = hp.read_map(f'{map_path}/map-f1z1.fits', dtype = int) # dtype makes sure it is integers

m = hp.ud_grade(m, 1024, power=-2) # this power preserves total number of objects

theo = np.loadtxt(f'{theo_path}/desy6-2404_f11z1f11z1.dat')

cl = hp.anafast(m)


print(cl.shape)
print('acabou')

plt.loglog(np.arange(3072), cl)
plt.loglog(np.arange(3072), m[np.arange(3072)])
plt.savefig('cl11.pdf')
