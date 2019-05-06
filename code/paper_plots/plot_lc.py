""" Collage of light curves for the Ic-BL SNe """

from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
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
dat = np.loadtxt(datadir + "/ztf.dat", dtype=str, delimiter=',')
names = dat[:,0]
nobj = len(names)
redshift = dat[:,1].astype(float)

# make a grid of say, 4 across
ncol = 4
nrow = int(nobj / ncol + 1)
fig,axarr = plt.subplots(nrow, ncol, figsize=(8,8), sharex=True, sharey=True)

# iterate across the grid
for ii,name in enumerate(names):
    ax = axarr.reshape(-1)[ii]
    name = names[ii]
    marshal.download_lightcurve(name)
    lc_dict = marshal.get_local_lightcurves(name)

    jd = lc_dict['jdobs'].values
    dm = Planck15.distmod(z=redshift[ii]).value
    mag = lc_dict['mag'].values
    emag = lc_dict['emag'].values
    filt = lc_dict['filter'].values

    choose = np.logical_and(filt == 'r', mag < 99)
    if sum(choose) > 0:
        ax.errorbar(
                (jd[choose]-jd[choose][0])/(1+redshift[ii]), mag[choose]-dm, 
                yerr=emag[choose], c='#140b34', fmt='s', ms=5)
    
    choose = np.logical_and(filt == 'g', mag < 99)
    if sum(choose) > 0:
        ax.errorbar(
                (jd[choose]-jd[choose][0])/(1+redshift[ii]), mag[choose]-dm, 
                yerr=emag[choose], c='#e55c30', fmt='o', ms=5)

    ax.text(0.95, 0.95, name.split("18")[1], transform=ax.transAxes, 
            fontsize=12, horizontalalignment='right',
            verticalalignment='top')
    ax.tick_params(labelsize=12, axis='both')

# Get rid of the extra panels
for ax in axarr.reshape(-1)[nobj:]:
    ax.set_visible(False)
    
fig.subplots_adjust(hspace=0, wspace=0)

fig.text(0.5, 0.04,
        "Rest-Frame Time Since Observed $r$-band Max [d]", fontsize=16,
        ha='center')

fig.text(0.04, 0.5, "Absolute Magnitude (AB)", 
        ha='center', va='center', fontsize=16, rotation='vertical')

ax.set_xlim(-20,100)
ax.invert_yaxis()
#plt.tight_layout()

plt.show()
#plt.savefig("lc_grid.eps", format='eps', dpi=1000)
