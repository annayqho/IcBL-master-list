""" Compile a list of Ic-BL SNe with z < 0.1 """

import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord,Distance
from astropy.cosmology import Planck15
from astropy.io import ascii


DATA_DIR = "/Users/annaho/Dropbox/Projects/Research/IcBL/data"


def todeg(ra, dec):
    """ convert XX:XX:XX to decimal degrees """
    radeg = []
    decdeg = []
    for ii,raval in enumerate(ra):
        hh = raval.split(":")[0]
        mm = raval.split(":")[1]
        ss = raval.split(":")[2]
        radegval = hh+"h"+mm+"m"+ss+"s"
        dd = dec[ii].split(":")[0]
        mm = dec[ii].split(":")[1]
        ss = dec[ii].split(":")[2]
        decdegval = dd+"d"+mm+"m"+ss+"s"
        c = SkyCoord(radegval, decdegval, frame='icrs')
        radeg.append(c.ra.deg)
        decdeg.append(c.dec.deg)
    return np.array(radeg), np.array(decdeg)


def ptf():
    """ the PTF/iPTF sample of 34 Ic-BL SNe
    I copied the table directly from the .tex file downloaded from the arXiv,
    then ran the following two commands
    %s/\\//g
    %s/ //g
    %s/\*//g
    %s/xx//g
    I also removed the commented-out lines
    """
    dat = Table.read(
            "%s/taddia2018.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    ra = dat['col2']
    dec = dat['col3']
    radeg, decdeg = todeg(ra, dec)
    z = dat['col5']
    return list(name), list(radeg), list(decdeg), list(z)


def ztf():
    """ The list of Ic-BL discovered in ZTF """
    dat = Table.read(
            "%s/ztf.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    ra = dat['col4']
    dec = dat['col5']
    radeg, decdeg = todeg(ra, dec)
    z = dat['col6']
    return list(name), list(radeg), list(decdeg), list(z)


def sdssII():
    """ the sample from SDSS-II (Taddia et al. 2015)
    I copied the table from the .tex file, and only kept the Ic-BL ones
    also got rid of the IIb that was commented out
    """
    dat = Table.read(
            "%s/taddia2015.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    ra_raw = dat['col2'].tolist()
    dec_raw = dat['col3'].tolist()
    for ii,ra in enumerate(ra_raw):
        ra_raw[ii] = ra.strip('$')
        dec_raw[ii] = dec_raw[ii].strip('$')
    ra = np.array(ra_raw)
    dec = np.array(dec_raw)
    z = np.array(dat['col5'])
    radeg, decdeg = todeg(ra, dec)
    return name, radeg, decdeg, z


def cano2013():
    """ The table of GRB/XRF-less Ic-BL SNe,
    from Cano et al. 2013 (since the ones with GRBs 
    I added positions to the table from the Open Supernova Catalog """
    dat = Table.read(
            "%s/cano2013.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    ra = dat['col3']
    dec = dat['col4']
    radeg,decdeg = todeg(ra,dec)
    z = dat['col5']
    return name, radeg, decdeg, z


def cano2016():
    """ The table of GRB-SNe from Cano et al. 2016 """
    dat = Table.read(
            "%s/cano2016.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    name = dat['col2'] # SN name, not GRB name
    ra = dat['col4']
    dec = dat['col5']
    radeg,decdeg = todeg(ra,dec)
    z = dat['col6']
    return name, radeg, decdeg, z


def lyman2016():
    """ The list of Ic-BL SNe from Lyman et al. 2016 """
    dat = Table.read(
            "%s/lyman2016.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    ra = dat['col3'] 
    dec = dat['col4'] 
    radeg,decdeg = todeg(ra,dec)
    temp = dat['col6'] 
    distmod = np.array(
            [val.split('pm')[0].strip('$') for val in temp]).astype(float)
    z = np.array([Distance(distmod=val).z for val in distmod])
    return name, radeg, decdeg, z


def prentice2016():
    """ The list of Ic-BL SNe from Prentice et al. 2016 """
    dat = Table.read(
            "%s/prentice2016.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    cl = dat['col2']
    ra = dat['col3']
    dec =dat['col4']
    radeg, decdeg = todeg(ra,dec)
    z = dat['col6']
    is_icbl = np.logical_or(cl=='Ic-BL', cl=='GRB-SN')
    return name[is_icbl], radeg[is_icbl], decdeg[is_icbl], z[is_icbl]


def modjaz2016():
    """ The list of Ic-BL SNe from Modjaz et al. 2016 
    Removed the one with a lower limit on redshift
    Actually removed all of them except SN2007bg, because
    that was the only new one. """
    dat = Table.read(
            "%s/modjaz.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    ra = dat['col2']
    dec =dat['col3']
    radeg, decdeg = todeg(ra,dec)
    z = dat['col4']
    return name, radeg, decdeg, z


def tns():
    """ This is the list on TNS as of 2018-12-29 """
    dat = Table.read(
            "%s/tns_search.csv" %DATA_DIR, 
            delimiter=';')
    name = dat['Name']
    ra = dat['RA']
    dec = dat['DEC']
    radeg, decdeg = todeg(ra,dec)
    z = dat['Redshift']
    return name, radeg, decdeg, z


def add(name, ra, dec, redshift, n, r, d, z):
    c = SkyCoord(ra, dec, unit='deg')
    cadd = SkyCoord(r, d, unit='deg')
    nadd = 0
    for ii,val in enumerate(cadd):
        dist = c.separation(val).arcsec
        nopos = False
        noname = False
        # Is the position in there already?
        if sum(dist <= 2) == 0: 
            nopos = True
        # Is the name in there already?
        if n[ii] not in name:
            noname = True
        if np.logical_and(nopos, noname):
            name.append(n[ii])
            ra.append(r[ii])
            dec.append(d[ii])
            redshift.append(z[ii])
            nadd += 1
        else:
            print("%s is a duplicate, not adding" %n[ii])
    print("added %s events" %str(nadd))
    return name, ra, dec, redshift


if __name__=="__main__":
    # Name, RA, Dec, Redshift 
    print("Adding the PTF/iPTF sample")
    name, ra, dec, redshift = ptf()
    print("added %s events" %len(name))

    # Add the new ZTF sample
    print("Adding the ZTF sample")
    n,r,d,z = ztf()
    name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    # Add the SDSS II sample
    print("Adding the SDSS II sample")
    n,r,d,z = sdssII()
    name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    # Add the Cano (2013) sample
    print("Adding the Cano (2013) sample")
    n,r,d,z = cano2013() 
    name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    # Add the Cano (2016) sample
    print("Adding the Cano (2016) sample")
    n,r,d,z = cano2016() 
    name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    # Add the Lyman (2016) sample
    print("Adding the Lyman (2016) sample")
    n,r,d,z = lyman2016() 
    name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    # Add the Prentice (2016) sample
    print("Adding the Prentice (2016) sample")
    n,r,d,z = prentice2016() 
    name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    # Add the Modjaz (2016) sample
    print("Adding the Modjaz (2016) sample")
    n,r,d,z = modjaz2016() 
    name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    # Add the most recent LLGRB (SN2017iuk)
    print("adding the most recent LLGRB")
    n = ['SN2017iuk']
    r = ['11:09:39.52']
    d = ['-12:35:18.34']
    radeg,decdeg = todeg(r,d)
    z = [0.037022]
    name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    # Add the list from TNS
    print("adding list from TNS")
    n,r,d,z = tns() 
    name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    name = np.array(name)
    redshift = np.array(redshift)
    ra = np.array(ra)
    dec = np.array(dec)
    print("TOTAL NUMBER IS")
    print(len(name))

    ascii.write(
            [name,ra,dec,redshift], 'all_icbl.html', 
            names=['Name', 'RA', 'Dec', 'z'], delimiter=',', overwrite=True,
            format='html')
