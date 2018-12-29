""" Check VLASS for Ic-BL SNe """

import sys
sys.path.append("/Users/annaho/Github/Query_VLASS")
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from vlass_search import search_vlass

# Get coordinates of all Ic-BL SNe
dat = ascii.read('all_icbl.html', format='html')
print(dat)
names = dat['Name']
ras = dat['RA']
decs = dat['Dec']
coords = SkyCoord(ras, decs, unit='deg')

# Check each one
for ii,name in enumerate(names):
    search_vlass(name, coords[ii])
