""" Collage of light curves for the Ic-BL SNe """

import numpy as np
from astropy.time import Time
from astropy.cosmology import Planck15
from ztfquery import query
from ztfquery import marshal

# connect to databases
m = marshal.MarshalAccess()
zquery = query.ZTFQuery()

# get the names of the Ic-BL SNe
datadir = "/Users/annaho/Dropbox/Projects/Research/IcBL/data"
names = np.loadtxt(datadir + "/ztf.dat", dtype=str)
nobj = len(names)

# get it working for one of them
name = names[0]
marshal.download_lightcurve(name)
lc_dict = marshal.get_local_lightcurves(name)
marshal.plot_lightcurve(lc_dict["marshal_lightcurve_ZTF18abcdef.csv"])
