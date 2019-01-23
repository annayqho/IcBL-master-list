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


def opensn():
    """ 
    Automatically grab all of the Ic-BL SNe from the open SN catalog """
    print("Connecting to the open SN catalog...")
    server = "https://api.sne.space/catalog"
    r = requests.get(server, params={'claimedtype': 'Ic BL', 'format': 'json'})
    dat = r.json()
    
    # Retrieve the data you want
    nsn = len(dat.keys())
    print("Found %s claimed Ic-BL SNe on the open SN catalog" %nsn)

    return dat


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
    dat = opensn()
    names = np.array(list(dat.keys()))
    nsn = len(names)
    ra = []
    dec = []
    for key,val in dat.items():
        if len(val['ra']) > 0:
            ra.append(val['ra'][0]['value'])
            dec.append(val['dec'][0]['value'])
    ra,dec = todeg(ra,dec)
    opensnpos = SkyCoord(ra, dec, unit='deg')

    # Question 1: are there any Ic-BL on TNS that are not on openSN?
    name, date, radeg, decdeg, z, ref = tns()
    name = np.array([val.replace(" ", "") for val in name])
    missing = np.setdiff1d(name,names)
    if len(missing) > 0:
        print("There are TNS Ic-BL SNe missing from OpenSN")
        print(missing)
    else:
        print("All TNS Ic-BL SNe are on OpenSN")

    # Question 2: are there any Ic-BL from other papers that are not on openSN?
    # Yes, a whole bunch from PTF and ZTF.
    name, date, radeg, decdeg, z, ref = ztf()
    name = np.array(name)
    print(np.setdiff1d(name,names))
    # compare positions, since some of these only have ZTF names...
    ptfpos = SkyCoord(radeg, decdeg, unit='deg')
    for ii,val in enumerate(ptfpos):
        if min(val.separation(opensnpos).arcsec) < 1:
            print("%s already in openSN" %name[ii])
        else:
            print("%s not in openSN" %name[ii])

    # # Name, Expl./Disc. Date, RA, Dec, Redshift, Reference
    # ascii.write(
    #         [names], 'all_icbl.html', names=['Name'], delimiter=',', 
    #         overwrite=True, format='html')
