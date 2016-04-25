import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from astropy.io import fits
import healpy as hp
import sys
import os

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

    return ra_arr,dec_arr,z0_arr,rsd_arr,type_arr

def mkmap(ra,dec,z0,rsd,zi,zf,nside,wrsd) :
    z=np.copy(z0)
    if wrsd :
        z+=rsd
    ind=np.where((z>=zi) & (z<zf))
    ra_here=ra[ind]
    dec_here=dec[ind]
    
    npix=hp.nside2npix(nside)
    mp_ng=np.histogram(hp.ang2pix(nside,np.pi*(90-dec_here)/180,ra_here*np.pi/180),bins=npix,range=[0,npix])[0]
    ndens=np.sum(mp_ng)/(4*np.pi)
    mp_delta=mp_ng/(np.mean(mp_ng)+0.0)-1
    
    return mp_delta,ndens

def run_colore(seed) :
    strout="prefix_out= dum\n"
    strout+="output_format= FITS\n"
    strout+="output_density= 0\n"
    strout+="pk_filename= outpk.dat\n"
    strout+="nz_filename= nz_test.txt\n"
    strout+="bias_filename= bias_test.txt\n"
    strout+="omega_M= 0.3\n"
    strout+="omega_L= 0.7\n"
    strout+="omega_B= 0.05\n"
    strout+="h= 0.7\n"
    strout+="w= -1.0\n"
    strout+="ns= 0.96\n"
    strout+="sigma_8= 0.824995\n"
    strout+="z_min= 0.05\n"
    strout+="z_max= 0.7\n"
    strout+="r_smooth= 2.\n"
    strout+="n_grid= 512\n"
    strout+="seed= %d\n"%seed
    f=open("param_dum.ini","w")
    f.write(strout)
    f.close()

    os.system('./CoLoRe param_dum.ini')
    ra_arr,dec_arr,z0_arr,rsd_arr,type_arr=read_fits("dum_0.fits")
    os.system('rm param_dum.ini dum_0.fits')

    return ra_arr,dec_arr,z0_arr,rsd_arr
              
if len(sys.argv)!= 6 :
    print "Usage: read_colore.py seed bins_file nside nsims prefix_out"
    exit(1)
seed=int(sys.argv[1])
fname_bin=sys.argv[2]
nside=int(sys.argv[3])
nsims=int(sys.argv[4])
prefix_out=sys.argv[5]
npix=hp.nside2npix(nside)
lmax=3*nside-1

zm_arr,dz_arr,sz_arr=np.loadtxt(fname_bin,unpack=True)
zi_arr=zm_arr-dz_arr
zf_arr=zm_arr+dz_arr
nbins=len(zi_arr)

cl_all=np.zeros([nsims,nbins,nbins,lmax+1])

for isim in np.arange(nsims) :
    print "Sim %d"%isim
    if os.path.isfile(prefix_out+"_s%d_cl_data.npy"%isim) : 
        cl_all[isim,:,:,:]=np.load(prefix_out+"_s%d_cl_data.npy"%isim)
    else :
        ra_arr,dec_arr,z0_arr,rsd_arr=run_colore(seed+isim)
        maps=np.zeros([nbins,npix])
        ndens=np.zeros(nbins)
        for ii in np.arange(nbins) :
            maps[ii,:],ndens[ii]=mkmap(ra_arr,dec_arr,z0_arr,rsd_arr,zi_arr[ii],zf_arr[ii],nside,True)
            print "Map %d"%ii+", shot noise = %lE"%(1./ndens[ii])
        for i1 in np.arange(nbins) :
            mp1=maps[i1,:]
            for i2 in np.arange(nbins-i1)+i1 :
                mp2=maps[i2,:]
                print "Cl %d-"%i1+"%d"%i2
                cl_all[isim,i1,i2,:]=hp.anafast(mp1,map2=mp2)
                if i1!=i2 :
                    cl_all[isim,i2,i1,:]=cl_all[isim,i1,i2,:]
            cl_all[isim,i1,i1,:]-=1./ndens[i1]
        np.save(prefix_out+"_s%d_cl_data.npy"%isim,cl_all[isim,:,:,:])

cl_mean=np.mean(cl_all,axis=0)
cl_rms=np.std(cl_all,axis=0)

larr=np.arange(lmax+1)
cols=['r','g','b','y']
data=np.loadtxt("outcl.dat",unpack=True)
larr_th=np.arange(len(data[0]))+2
cl_th=np.zeros([nbins,nbins,len(data[0])])
icl=0
for i1 in np.arange(nbins) :
    for i2 in np.arange(nbins-i1)+i1 :
        cl_th[i1,i2,:]=data[icl+1]*2*np.pi/(larr_th*(larr_th+1))
        if i1!=i2 :
            cl_th[i2,i1,:]=cl_th[i1,i2,:]
        icl+=1

for i1 in np.arange(nbins) :
    label="Bin %d"%(i1+1)+", $z\\in[%.2lf,"%(zi_arr[i1])+"%.2lf)$"%(zf_arr[i1])
    plt.errorbar(larr,cl_mean[i1,i1],yerr=cl_rms[i1,i1]/np.sqrt(nsims),fmt=cols[i1]+'.',label=label)
    plt.plot(larr_th,cl_th[i1,i1],cols[i1]+'-')
plt.errorbar(larr,-np.ones_like(larr),yerr=cl_rms[0,0],fmt='k.',label='CoLoRe')
plt.plot(larr,-np.ones_like(larr),'k-',label='Linear theory')
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.ylim([1E-6,2E-3])
plt.xlim([2,lmax])
plt.xlabel('$\\ell$',fontsize=16)
plt.ylabel('$C_\\ell$',fontsize=16)
plt.legend(loc='upper right',frameon=False,labelspacing=0.1)
plt.show()
