""" Check VLASS for Ic-BL SNe """

import numpy as np
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
ref = {}
for ii,name in enumerate(names):
    ref[name] = coords[ii]

# Check each one
for name in names:
    if name != 'iPTF15dqg':
        search_vlass(name, ref[name])
