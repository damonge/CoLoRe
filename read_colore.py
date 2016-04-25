import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from astropy.io import fits
import healpy as hp
import sys

def read_ascii(fname) :
    ra_arr,dec_arr,z0_arr,rsd_arr,type_arr=np.loadtxt(fname,unpack=True)

    return ra_arr,dec_arr,z0_arr+rsd_arr,type_arr

def read_hdf5(fname) :
    ff=h5.File(fname,"r")
    
    ra_arr=ff['/sources']['RA']
    dec_arr=ff['/sources']['DEC']
    z0_arr=ff['/sources']['Z_COSMO']
    rsd_arr=ff['/sources']['DZ_RSD']
    type_arr=ff['/sources']['TYPE']

    return ra_arr,dec_arr,z0_arr+rsd_arr,type_arr

def read_fits(fname) :
    data=(fits.open(fname)[1]).data
    
    ra_arr=data['RA']
    dec_arr=data['DEC']
    z0_arr=data['Z_COSMO']
    rsd_arr=data['DZ_RSD']
    type_arr=data['TYPE']

    return ra_arr,dec_arr,z0_arr+rsd_arr,type_arr

if len(sys.argv)!= 3 :
    print "Usage: read_colore.py file_name file_format ('ASCII', 'FITS' or 'HDF5')"
    exit(1)
fname=sys.argv[1]
fmt=sys.argv[2]

if fmt=='ASCII' :
    ra_arr,dec_arr,z_arr,type_arr=read_ascii(fname)
elif fmt=='FITS' :
    ra_arr,dec_arr,z_arr,type_arr=read_fits(fname)
elif fmt=='HDF5' :
    ra_arr,dec_arr,z_arr,type_arr=read_hdf5(fname)

nside=64
npix=hp.nside2npix(nside)
mp=np.histogram(hp.ang2pix(nside,np.pi*(90-dec_arr[np.where(type_arr==0)])/180,np.pi*ra_arr[np.where(type_arr==0)]/180),bins=npix,range=[0,npix])[0]
plt.figure();
plt.hist(z_arr[np.where(type_arr==0)],bins=100,histtype='step');
plt.hist(z_arr[np.where(type_arr==1)],bins=100,histtype='step');
plt.xlabel('$z$',fontsize=16); plt.ylabel('$N(z)$',fontsize=16);
hp.mollview(mp,title='$N_g(\\hat{\\bf n})$');
plt.show()
