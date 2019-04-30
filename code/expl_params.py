""" Play around with explosion parameters """

import numpy as np

dat = np.loadtxt(
    "/Users/annaho/Dropbox/Projects/Research/IcBL/data/expl_params.dat",
    delimiter='&', dtype=str)

names = dat[:,0]
nobj = len(names)

mej_raw = dat[:,1]
mej = np.zeros(nobj)
emej = np.zeros(nobj)

ek_raw = dat[:,2]
ek = np.zeros(nobj)
eek = np.zeros(nobj)

mni_raw = dat[:,3]
mni = np.zeros(nobj)
emni = np.zeros(nobj)

mixing = dat[:,4]

for ii,val in enumerate(mej_raw):
    mej[ii] = float(mej_raw[ii].split('(')[0])
    emej[ii] = float(mej_raw[ii].split('(')[1].split(')')[0])

for ii,val in enumerate(ek_raw):
    ek[ii] = float(ek_raw[ii].split('(')[0])
    eek[ii] = float(ek_raw[ii].split('(')[1].split(')')[0])

for ii,val in enumerate(mni_raw):
    mni[ii] = float(mni_raw[ii].split('(')[0])
    emni[ii] = float(mni_raw[ii].split('(')[1].split(')')[0])


