import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from astropy.io import fits
import healpy as hp
import sys
import os

def read_ascii(fname) :
    ifile=0
    
    type_arr,ra_arr,dec_arr,z0_arr,rsd_arr,e1_arr,e2_arr=np.loadtxt(fname+"_%d.txt"%ifile,unpack=True)

    ifile=1
    while os.path.isfile(fname+"_%d.h5"%ifile) :
        typ,ra,dec,z0,rsd,e1,e2=np.loadtxt(fname+"_%d.txt"%ifile,unpack=True)
        type_arr=np.concatenate((type_arr,typ_arr))
        ra_arr=np.concatenate((ra_arr,ra_arr))
        dec_arr=np.concatenate((dec_arr,dec_arr))
        z0_arr=np.concatenate((z0_arr,z0_arr))
        rsd_arr=np.concatenate((rsd_arr,rsd_arr))
        e1_arr=np.concatenate((e1_arr,e1_arr))
        e2_arr=np.concatenate((e2_arr,e2_arr))
        ifile+=1

    return ra_arr,dec_arr,z0_arr,rsd_arr,type_arr,e1_arr,e2_arr

def read_hdf5(fname,ipop) :
    ifile=0
    
    ff=h5.File(fname+"_%d.h5"%ifile,"r")
    tabname='/sources%d'%ipop

    ra_arr=ff[tabname]['RA']
    dec_arr=ff[tabname]['DEC']
    z0_arr=ff[tabname]['Z_COSMO']
    rsd_arr=ff[tabname]['DZ_RSD']
    e1_arr=ff[tabname]['E1']
    e2_arr=ff[tabname]['E2']

    ifile=1
    while os.path.isfile(fname+"_%d.h5"%ifile) :
        ff=h5.File(fname+"_%d.h5"%ifile,"r")
        tabname='/sources%d'%ipop
        
        ra_arr=np.concatenate((ra_arr,ff[tabname]['RA']))
        dec_arr=np.concatenate((dec_arr,ff[tabname]['DEC']))
        z0_arr=np.concatenate((z0_arr,ff[tabname]['Z_COSMO']))
        rsd_arr=np.concatenate((rsd_arr,ff[tabname]['DZ_RSD']))
        e1_arr=np.concatenate((e1_arr,ff[tabname]['E1']))
        e2_arr=np.concatenate((e2_arr,ff[tabname]['E2']))
        ifile+=1

    type_arr=ipop*np.ones(len(ra_arr))

    return ra_arr,dec_arr,z0_arr,rsd_arr,type_arr,e1_arr,e2_arr

def read_fits(fname) :
    ifile=0
    
    data=(fits.open(fname+"_%d.fits"%ifile)[1]).data
    
    ra_arr=data['RA']
    dec_arr=data['DEC']
    z0_arr=data['Z_COSMO']
    rsd_arr=data['DZ_RSD']
    e1_arr=data['E1']
    e2_arr=data['E2']
    type_arr=data['TYPE']

    ifile=1
    while os.path.isfile(fname+"_%d.fits"%ifile) :
        data=(fits.open(fname+"_%d.fits"%ifile)[1]).data

        ra_arr=np.concatenate((ra_arr,data['RA']))
        dec_arr=np.concatenate((dec_arr,data['DEC']))
        z0_arr=np.concatenate((z0_arr,data['Z_COSMO']))
        rsd_arr=np.concatenate((rsd_arr,data['DZ_RSD']))
        e1_arr=np.concatenate((e1_arr,data['E1']))
        e2_arr=np.concatenate((e2_arr,data['E2']))
        type_arr=np.concatenate((type_arr,data['TYPE']))
        ifile+=1

    return ra_arr,dec_arr,z0_arr,rsd_arr,type_arr,e1_arr,e2_arr

if len(sys.argv)!= 3 :
    print "Usage: read_colore.py file_name file_format ('ASCII', 'FITS' or 'HDF5')"
    exit(1)
fname=sys.argv[1]
fmt=sys.argv[2]

if fmt=='ASCII' :
    ra_arr,dec_arr,z_arr,rsd_arr,type_arr,e1_arr,e2_arr=read_ascii(fname)
elif fmt=='FITS' :
    ra_arr,dec_arr,z_arr,rsd_arr,type_arr,e1_arr,e2_arr=read_fits(fname)
elif fmt=='HDF5' :
    ra_arr,dec_arr,z_arr,rsd_arr,type_arr,e1_arr,e2_arr=read_hdf5(fname,0)

nside=64
npix=hp.nside2npix(nside)
mp=np.histogram(hp.ang2pix(nside,np.pi*(90-dec_arr[np.where(type_arr==0)])/180,np.pi*ra_arr[np.where(type_arr==0)]/180),bins=npix,range=[0,npix])[0]
mpr=np.histogram(hp.ang2pix(nside,np.pi*(90-dec_arr[np.where(type_arr==0)])/180,np.pi*ra_arr[np.where(type_arr==0)]/180),bins=npix,range=[0,npix],weights=rsd_arr)[0]
mpe1=np.histogram(hp.ang2pix(nside,np.pi*(90-dec_arr[np.where(type_arr==0)])/180,np.pi*ra_arr[np.where(type_arr==0)]/180),bins=npix,range=[0,npix],weights=e1_arr)[0]
mpe2=np.histogram(hp.ang2pix(nside,np.pi*(90-dec_arr[np.where(type_arr==0)])/180,np.pi*ra_arr[np.where(type_arr==0)]/180),bins=npix,range=[0,npix],weights=e2_arr)[0]
plt.figure();
plt.hist(z_arr[np.where(type_arr==0)],bins=100,histtype='step');
plt.hist(z_arr[np.where(type_arr==1)],bins=100,histtype='step');
plt.xlabel('$z$',fontsize=16); plt.ylabel('$N(z)$',fontsize=16);
hp.mollview(mp,title='$N_g(\\hat{\\bf n})$')
hp.mollview(mpr/mp,title='$v_r$')
hp.mollview(mpe1/mp,title='$e_1$')
hp.mollview(mpe2/mp,title='$e_2$')
plt.show()
