""" 
Print a table of the Ic-BL sample from ZTF
"""

import numpy as np
from math import floor, log10
from astropy.time import Time
from astropy.cosmology import Planck15
from ztfquery import query
from ztfquery import marshal


def round_sig(x, sig=2):
    if x < 0:
        return -round(-x, sig-int(floor(log10(-x)))-1)
    return round(x, sig-int(floor(log10(x)))-1)


def ndec(num):
    dec = str(num).split('.')[-1]
    return len(dec)

# connect to databases
m = marshal.MarshalAccess()
zquery = query.ZTFQuery()
m.load_target_sources()

# initialize LaTeX table
headings = np.array(
        ['ZTF ID', 'RA (hms)', 'Dec (dms)', 'Disc (JD)', 
         '$z$', 'Ref.'])
label = "sample"
caption = "Sample of Ic-BL SNe from ZTF in the period May 2018 -- May 2019"

# Print the table headers
ncol = len(headings)
colstr = ""
colstr += 'l'
for col in np.arange(ncol-1): colstr+="r"
print(colstr)

colheadstr = ""
for col in np.arange(ncol-1):
    colheadstr += "\colhead{%s} & " %headings[col]
colheadstr += "\colhead{%s}" %headings[-1]

rowstr = ""
for col in np.arange(ncol-1):
    rowstr += "%s & "
rowstr += "%s \\\ \n"

outputf = open("table_%s.txt" %label, "w")
outputf.write("\\startlongtable \n")
outputf.write("\\begin{deluxetable}{%s} \n" %colstr)
outputf.write("\\tablecaption{%s\label{tab:%s}} \n" %(caption,label))
outputf.write("\\tablewidth{0pt} \n")
outputf.write("\\tablehead{ %s } \n" %colheadstr)
#outputf.write("\\rotate \n")
outputf.write("\\tabletypesize{\scriptsize} \n")
outputf.write("\startdata \n")

datadir = "/Users/annaho/Dropbox/Projects/Research/IcBL/data"
names = np.loadtxt(datadir + "/ztf.dat", dtype=str)
nobj = len(names)

# extract the redshifts
redshift = m.get_target_redshift(names).values.astype(float)

for ii in np.arange(nobj):
    name = names[ii]
    row = rowstr %(name, names[ii], names[ii], names[ii], redshift[ii], names[ii])
    outputf.write(row)

outputf.write("\enddata \n")
outputf.write("\end{deluxetable} \n")
