""" Redshift distribution of all Ic-BL """

from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
from astropy.table import Table

dat = Table.read("all_icbl.html", format='html')
z = dat['z']

fig = plt.figure(figsize=(5,3))

plt.hist(z, histtype='step', range=(0,1.01), color='black', lw=2)
plt.title("$z$ distribution for all Ic-BL SNe", fontsize=14)
plt.xlabel("Redshift", fontsize=14)
plt.ylabel("Count", fontsize=14)
plt.tick_params(axis='both', labelsize=14)

plt.tight_layout()

#plt.show()
plt.savefig("icbl_z_dist.png")
