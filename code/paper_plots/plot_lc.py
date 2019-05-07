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
import extinction


# Global variables
datadir = "/Users/annaho/Dropbox/Projects/Research/IcBL/data"


def plot_98bw(ax):
    # Plot the light curve of SN1998bw
    dm = Planck15.distmod(0.0085).value
    lc = np.loadtxt(datadir + "/sn1998bw.dat", delimiter=';', dtype=str)
    jd = lc[:,0].astype(float)
    rmag_raw = lc[:,7]
    choose = rmag_raw != '     '
    rmag = rmag_raw[choose]
    zp = jd[choose][np.argmin(rmag)]
    ax.plot(
            (jd[choose]-zp)/(1.0085),
            rmag.astype(float)-dm, lw=0.5, c='black')
    vmag_raw = lc[:,5]
    choose = vmag_raw != '     '
    vmag = vmag_raw[choose] 
    ax.plot(
            (jd[choose]-zp)/(1.0085), 
            vmag.astype(float)-dm, lw=0.5, c='#e55c30', ls='--')

# connect to databases
m = marshal.MarshalAccess()
zquery = query.ZTFQuery()
print("Connected")

# download metadata for all sources
m.load_target_sources()
print("Downloaded metadata")

# get the names of the Ic-BL SNe
dat = np.loadtxt(datadir + "/ztf.dat", dtype=str, delimiter=',')
names = dat[:,0]

# get redshifts
redshift = m.get_target_redshift(names).values.astype(float)

# make a grid of say, 4 across
nobj = len(names)
ncol = 3
nrow = int(nobj / ncol + 1)
fig,axarr = plt.subplots(nrow, ncol, figsize=(7,9), sharex=True, sharey=True)

# iterate across the grid
for ii,name in enumerate(names):
    ax = axarr.reshape(-1)[ii]
    plot_98bw(ax)

    name = names[ii]

    # Get the light curve
    marshal.download_lightcurve(name)
    lc_dict = marshal.get_local_lightcurves(name)

    # Get the epochs of spectra
    spec_fnames = list(marshal.download_spectra(name, dirout=None).keys())
    specdates_raw = np.array([f.split('_')[1] for f in spec_fnames])
    specdates = np.array(
            ['-'.join([d[0:4],d[4:6],d[6:]]) for d in specdates_raw])
    specjd = Time(specdates).jd

    # Plot the light curve
    jd = lc_dict['jdobs'].values
    dm = Planck15.distmod(z=redshift[ii]).value
    mag = lc_dict['mag'].values
    emag = lc_dict['emag'].values
    filt = lc_dict['filter'].values

    choose = np.logical_and(filt == 'r', mag < 99)
    zp = jd[choose][np.argmin(mag[choose])]

    if sum(choose) > 0:
        ax.errorbar(
                (jd[choose]-zp)/(1+redshift[ii]), mag[choose]-dm, 
                yerr=emag[choose], c='#140b34', fmt='s', ms=5)
    
    choose = np.logical_and(filt == 'g', mag < 99)
    if sum(choose) > 0:
        ax.errorbar(
                (jd[choose]-zp)/(1+redshift[ii]), mag[choose]-dm, 
                yerr=emag[choose], c='#e55c30', fmt='o', ms=5)

    # Plot the epochs of spectra
    for t in specjd:
        ax.axvline(x=t-zp, lw=1.0, c='lightgrey', ls='-')

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

ax.set_xlim(-30,100)
ax.set_ylim(-15,-20.5)
#ax.invert_yaxis()
#plt.tight_layout()

#plt.show()
plt.savefig("lc_grid.eps", format='eps', dpi=1000)
