import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
import os.path

nside=128
npix=hp.nside2npix(nside)

#Analyze sources
ifile=0
nmap=np.zeros(hp.nside2npix(nside))
e1map=np.zeros(hp.nside2npix(nside))
e2map=np.zeros(hp.nside2npix(nside))
while os.path.isfile('test/out_srcs_%d.fits'%ifile) :
    hdulist = fits.open('test/out_srcs_%d.fits'%ifile)
    tbdata = hdulist[1].data
    
    pix=hp.ang2pix(nside,(90-tbdata['DEC'])*np.pi/180,tbdata['RA']*np.pi/180)
    n=np.bincount(pix,minlength=npix)
    e1=np.bincount(pix,minlength=npix,weights=tbdata['E1'])
    e2=np.bincount(pix,minlength=npix,weights=tbdata['E2'])
    nmap+=n
    e1map+=e1
    e2map+=e2
    ifile+=1

ndens=(np.sum(nmap)+0.0)/(4*np.pi)
mp_e1=e1map/nmap; mp_e1[nmap<=0]=0
mp_e2=e2map/nmap; mp_e2[nmap<=0]=0
mp_d=(nmap+0.0)/np.mean(nmap+0.0)-1
mp_db,mp_E,mp_B=hp.alm2map(hp.map2alm(np.array([mp_d,mp_e1,mp_e2]),pol=True),pol=False,nside=nside)

lt,cls_dd=np.loadtxt('test/outlj_cl_dd.txt',unpack=True);
lt,clt_kd=np.loadtxt('test/outlj_cl_dc.txt',unpack=True);
lt,clt_kk=np.loadtxt('test/outlj_cl_cc.txt',unpack=True);
lt,clt_dl=np.loadtxt('test/outlj_cl_d1l2.txt',unpack=True);
lt,clt_ll=np.loadtxt('test/outlj_cl_ll.txt',unpack=True);
cln_dd=np.ones_like(lt)/ndens
clt_dd=cls_dd+cln_dd
cld_dd,cld_ee,cld_bb,cld_de,cld_eb,cld_db=hp.anafast(np.array([mp_d,mp_e1,mp_e2]),pol=True);
ld=np.arange(len(cld_dd));

#Analyze kappa
mp_k=hp.read_map("test/out_kappa_z000.fits")
cld_kk=hp.anafast(mp_k); ld=np.arange(len(cld_kk))
cld_kd=hp.anafast(mp_k,map2=mp_d)

#Plots
hp.mollview(mp_d);
hp.mollview(mp_E);
hp.mollview(mp_B);
hp.mollview(mp_e1);
hp.mollview(mp_e2);
hp.mollview(mp_k)
plt.figure()
plt.hist(mp_e1,bins=100,histtype='step');
plt.hist(mp_e2,bins=100,histtype='step');
plt.figure()
plt.plot(ld,cld_dd,'r-',label='$\\delta_g\\times\\delta_g$',lw=2)
plt.plot(ld,cld_ee,'y-',label='$\\gamma^E_g\\times\\gamma^E_g$',lw=2)
plt.plot(ld,cld_de,'c-',label='$\\gamma^E_g\\times\\delta_g$',lw=2)
plt.plot(ld,cld_bb,'m',linestyle='solid' ,label='$\\gamma^B_g\\times\\gamma^B_g$',lw=2)
plt.plot(ld,cld_eb,'m',linestyle='dashed',label='$\\gamma^E_g\\times\\gamma^B_g$',lw=2)
plt.plot(ld,cld_db,'m',linestyle='dotted',label='$\\gamma^B_g\\times\\delta_g$',lw=2)
plt.plot(ld,cld_kd,'g-',label='$\\kappa-\\delta_g$',lw=2)
plt.plot(ld,cld_kk,'b-',label='$\\kappa-\\kappa$',lw=2)
plt.plot(lt,clt_dd,'k-',lw=2)
plt.plot(lt,2*clt_dl,'k-',lw=2)
plt.plot(lt,4*clt_ll,'k-',lw=2)
plt.plot(lt,clt_kd,'k-',lw=2)
plt.plot(lt,clt_kk,'k-',lw=2)
plt.loglog()
plt.xlabel('$\\ell$',fontsize=16)
plt.ylabel('$C_\\ell$',fontsize=16)
plt.xlim([2,192])
plt.legend(loc='lower left',frameon=False,labelspacing=0.1,ncol=2)
plt.savefig("test/cls.png",bbox_inches='tight')
plt.show()
