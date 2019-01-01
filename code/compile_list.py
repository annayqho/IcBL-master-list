""" Compile a list of Ic-BL SNe """

import numpy as np
import requests
from astropy.table import Table
from astropy.time import Time
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


def openSN():
    """ 
    Automatically grab all of the Ic-BL SNe from the open SN catalog """
    print("Connecting to the open SN catalog...")
    server = "https://api.sne.space/catalog"
    r = requests.get(server, params={'claimedtype': 'Ic BL', 'format': 'csv'})
    alldat = r.text.split('\n')
    
    # Header
    header = np.array(alldat[0].split(','))

    # Data
    dat = alldat[1:]
    # According to the formatting, you want to group things that live together
    # in double quotation marks. So, the real split between items is ",", not ,
    for ii,row in enumerate(dat):
        dat[ii] = np.array(dat[ii].split('","'))
    dat = np.array(dat)

    # Retrieve the data you want
    nsn = dat.shape[0]
    print("Found %s claimed Ic-BL SNe on the open SN catalog" %nsn)

    name = dat[:,np.where(header=='Name')[0][0]]
    ra = dat[:,np.where(header=='RA')[0][0]]
    dec = dat[:,np.where(header=='DEC')[0][0]]
    radeg, decdeg = todeg(ra,dec)
    z = dat[:,np.where(header=='Redshift')[0][0]]
    date = dat[:,np.where(header=='Discovery Date (UT)')[0][0]]
    ref = ['TNS'] * nsn

    return name, date, radeg, decdeg, z, ref


def tns():
    """ Run this to automatically grab all of the Ic-BL SNe from TNS """
    print("Connecting to TNS server...")
    server = "https://wis-tns.weizmann.ac.il/search"
    r = requests.get(server, params={'objtype': 7, 'format': 'csv'})
    alldat = r.text.split('\n')
    
    # Header
    header = np.array(alldat[0].split('","'))

    # Data
    dat = alldat[1:]
    # According to the formatting, you want to group things that live together
    # in double quotation marks. So, the real split between items is ",", not ,
    for ii,row in enumerate(dat):
        dat[ii] = np.array(dat[ii].split('","'))
    dat = np.array(dat)

    # Retrieve the data you want
    nsn = dat.shape[0]
    print("Found %s Ic-BL SNe on TNS" %nsn)

    name = dat[:,np.where(header=='Name')[0][0]]
    ra = dat[:,np.where(header=='RA')[0][0]]
    dec = dat[:,np.where(header=='DEC')[0][0]]
    radeg, decdeg = todeg(ra,dec)
    z = dat[:,np.where(header=='Redshift')[0][0]]
    date = dat[:,np.where(header=='Discovery Date (UT)')[0][0]]
    ref = ['TNS'] * nsn

    return name, date, radeg, decdeg, z, ref



def ptf():
    """ the PTF/iPTF sample of 34 Ic-BL SNe
    I copied the table directly from the .tex file downloaded from the arXiv,
    then ran the following two commands
    %s/\\//g
    %s/ //g
    %s/\*//g
    %s/xx//g
    I also removed the commented-out lines
    In this paper, they give estimated explosion epochs (with a typical
    uncertainty of 2 days) for all of the SNe observed before
    and after r maximum brightness.
    A lot of them don't have an estimated explosion epoch, though.
    So what I should do is use the estimate for the ones that have it,
    and for the ones that don't have it, just report discovery date
    as I found it on the marshal.
    """
    # Discovery dates on the Marshal, for the ones that aren't in Table 2
    # 27 out of 34 leaves 7
    disc = {}
    disc['PTF09sk'] = 2455002.74571
    disc['PTF10cs'] = 2455203.74537
    disc['PTF12grr'] = 2456117.84878
    disc['iPTF14bfu'] = Time('2014-06-06T03:11:51.86').jd
    disc['iPTF15dld'] = 2457318.82184
    disc['iPTF16coi'] = 2457625.72566
    disc['iPTF17axg'] = 2457784.97286

    dat = Table.read(
            "%s/taddia2018.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    # file with explosion epochs
    dat_expl = Table.read(
            "%s/taddia2018_t2.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    name_expl = dat_expl['col1']
    texpl = dat_expl['col8']

    name = dat['col1']
    texpl = []
    for n in name:
        try: 
            ind = np.where(name_expl==n)[0][0]
            texpl.append(texpl[ind])
        except:
            texpl.append(disc[n])
    ra = dat['col2']
    dec = dat['col3']
    radeg, decdeg = todeg(ra, dec)
    z = dat['col5']

    ref = ['T18']*len(name)

    return list(name), texpl, list(radeg), list(decdeg), list(z), ref


def ztf():
    """ The list of Ic-BL discovered in ZTF """
    dat = Table.read(
            "%s/ztf.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    name = dat['col1']
    date = dat['col3']
    ra = dat['col5']
    dec = dat['col6']
    radeg, decdeg = todeg(ra, dec)
    z = dat['col7']
    ref = ['ZTF']*len(name)
    return list(name), list(date), list(radeg), list(decdeg), list(z), ref


def sdssII():
    """ the sample from SDSS-II (Taddia et al. 2015)
    I copied the table from the .tex file, and only kept the Ic-BL ones
    also got rid of the IIb that was commented out
    """
    dat = Table.read(
            "%s/taddia2015.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    # Discovery dates from TNS
    disc = {}
    disc['SN2005fk'] = Time('2005-09-12T00:00:00').jd
    disc['SN2005kr'] = Time('2005-11-03T00:00:00').jd
    # from open SN catalog
    disc['SN2005ks'] = Time('2005-11-04T00:00:00').jd
    # from Table 3 in Taddia 2015...not sure what t_0 is though
    disc['SN14475'] = Time(54008.68, format='mjd').jd
    # actual texpl from Taddia 2015
    disc['SN2006nx'] = Time(54038.89, format='mjd').jd

    name = dat['col1']
    ra_raw = dat['col2'].tolist()
    dec_raw = dat['col3'].tolist()
    texpl = []
    ref = []
    for ii,ra in enumerate(ra_raw):
        ra_raw[ii] = ra.strip('$')
        dec_raw[ii] = dec_raw[ii].strip('$')
        texpl.append(disc[name[ii]])
        ref.append('T15')
    ra = np.array(ra_raw)
    dec = np.array(dec_raw)
    z = np.array(dat['col5'])
    radeg, decdeg = todeg(ra, dec)
    return name, texpl, radeg, decdeg, z, ref


def cano2013():
    """ The table of GRB/XRF-less Ic-BL SNe,
    from Cano et al. 2013 (since the ones with GRBs 
    I added positions from the Open Supernova Catalog """
    dat = Table.read(
            "%s/cano2013.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    disc = {}
    disc['SN1997ef'] = Time('1997-11-25T00:00:00').jd
    disc['SN2002ap'] = Time('2002-01-29T00:00:00').jd
    disc['SN2003jd'] = Time('2003-10-25T00:00:00').jd
    disc['SN2005kz'] = Time('2005-12-01T00:00:00').jd
    disc['SN2007D'] = Time('2007-01-09T00:00:00').jd
    disc['SN2007ru'] = Time('2007-11-27T00:00:00').jd
    disc['SN2009bb'] = Time('2009-03-21T00:00:00').jd
    disc['SN2010ah'] = Time('2010-02-23T00:00:00').jd
    disc['SN2010ay'] = Time('2010-03-17T00:00:00').jd
    name = dat['col1']
    ra = dat['col3']
    dec = dat['col4']
    radeg,decdeg = todeg(ra,dec)
    z = dat['col5']
    date = []
    ref = []
    for n in name:
        date.append(disc[n])
        ref.append('C13')
    return name, date, radeg, decdeg, z, ref


def cano2016():
    """ The table of GRB-SNe from Cano et al. 2016 """
    dat = Table.read(
            "%s/cano2016.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    grbname = dat['col1']
    date = []
    # turn GRB name into a date
    for n in grbname:
        yy = n[0:2]
        mm = n[2:4]
        dd = n[4:6]
        if yy[0] == '9':
            tstr = '19%s-%s-%sT00:00:00' %(yy,mm,dd)
            date.append(Time(tstr, format='isot').jd)
        else:
            tstr = '20%s-%s-%sT00:00:00' %(yy,mm,dd)
            date.append(Time(tstr, format='isot').jd)
    name = dat['col2'] # SN name
    ra = dat['col4']
    dec = dat['col5']
    radeg,decdeg = todeg(ra,dec)
    z = dat['col6']
    ref = ['C16'] * len(name)
    return name, date, radeg, decdeg, z, ref


def lyman2016():
    """ The list of Ic-BL SNe from Lyman et al. 2016 """
    dat = Table.read(
            "%s/lyman2016.dat" %DATA_DIR, 
            delimiter='&', format='ascii.fast_no_header')
    # from open SN
    disc = {}
    disc['SN1998bw'] = Time('1998-04-28T00:00:00').jd
    disc['SN2002ap'] = Time('2002-01-29T00:00:00').jd
    disc['SN2003jd'] = Time('2003-10-25T00:00:00').jd
    disc['SN2005kz'] = Time('2005-12-01T00:00:00').jd
    disc['SN2006aj'] = Time('2006-02-18T00:00:00').jd
    disc['SN2007ru'] = Time('2007-11-27T00:00:00').jd
    disc['SN2009bb'] = Time('2009-03-21T00:00:00').jd
    disc['SN2010bh'] = Time('2010-03-16T00:00:00').jd
    date = []
    name = dat['col1']
    for n in name:
       date.append(disc[n]) 
    ra = dat['col3'] 
    dec = dat['col4'] 
    radeg,decdeg = todeg(ra,dec)
    temp = dat['col6'] 
    distmod = np.array(
            [val.split('pm')[0].strip('$') for val in temp]).astype(float)
    z = np.array([Distance(distmod=val).z for val in distmod])
    ref = ['L16']*len(name)
    return name, date, radeg, decdeg, z, ref


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



def add(name, disc, ra, dec, redshift, ref, n, di, r, d, z, re):
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
            disc.append(di[ii])
            ra.append(r[ii])
            dec.append(d[ii])
            redshift.append(z[ii])
            ref.append(re[ii])
            nadd += 1
        else:
            print("%s is a duplicate, not adding" %n[ii])
    print("added %s events" %str(nadd))
    return name, disc, ra, dec, redshift, ref


if __name__=="__main__":
    # Name, Approximate Expl. Date, RA, Dec, Redshift, Reference

    print("Adding the PTF/iPTF sample")
    name, disc, ra, dec, redshift, ref = ptf()
    print("added %s events" %len(name))

    # Add the new ZTF sample
    print("Adding the ZTF sample")
    n,di,r,d,z,re = ztf()
    name, disc, ra, dec, redshift, ref = add(
            name, disc, ra, dec, redshift, ref, n, di, r, d, z, re)

    # Add the SDSS II sample
    print("Adding the SDSS II sample")
    n,di,r,d,z,re = sdssII()
    name, disc, ra, dec, redshift, ref = add(
            name, disc, ra, dec, redshift, ref, n, di, r, d, z, re)

    # Add the Cano (2013) sample
    print("Adding the Cano (2013) sample")
    n,di,r,d,z,re = cano2013() 
    name, disc, ra, dec, redshift, ref = add(
            name, disc, ra, dec, redshift, ref, n, di, r, d, z, re)

    # Add the Cano (2016) sample
    print("Adding the Cano (2016) sample")
    n,di,r,d,z,re = cano2016() 
    name, disc, ra, dec, redshift, ref = add(
            name, disc, ra, dec, redshift, ref, n, di, r, d, z, re)

    # Add the Lyman (2016) sample
    print("Adding the Lyman (2016) sample")
    n,di,r,d,z,re = lyman2016() 
    name, disc, ra, dec, redshift, ref = add(
            name, disc, ra, dec, redshift, ref, n, di, r, d, z, re)

    # Add the Prentice (2016) sample
    print("Adding the Prentice (2016) sample")
    n,r,d,z = prentice2016() 
    print(n)
    # name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    # Add the Modjaz (2016) sample
    # print("Adding the Modjaz (2016) sample")
    # n,r,d,z = modjaz2016() 
    # name, ra, dec, redshift = add(name, ra, dec, redshift, n, r, d, z)

    # Add the list from TNS
    print("adding list from TNS")
    n,di,r,d,z,re = tns() 
    name, disc, ra, dec, redshift, ref = add(
            name, disc, ra, dec, redshift, ref, n, di, r, d, z, re)

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
