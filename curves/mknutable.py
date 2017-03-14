import numpy as np
import py_cosmo_mad as csm
from scipy.optimize import brentq
from scipy.interpolate import interp1d

pcs=csm.PcsPar()
pcs.background_set(0.3,0.7,0.05,-1.,0.,0.7,2.7255)

def get_nu_next(nu0,dr) :
    a0=nu0/1420.
    r0=pcs.radial_comoving_distance(a0)
    def f_min(a) :
        return pcs.radial_comoving_distance(a)+dr-r0
    a_next=brentq(f_min,0.1,1.)
    
    return a_next*1420.

nu0_arr=[]
nuf_arr=[]

nu0=740.
while nu0<=1049. :
    nu0_arr.append(nu0)
    nu0=get_nu_next(nu0,5.)
    nuf_arr.append(nu0)
nu0_arr=np.array(nu0_arr)
nuf_arr=np.array(nuf_arr)

np.savetxt("nus.txt",np.transpose([nu0_arr,nuf_arr]))

zDESI_in,nzDESI_in=np.loadtxt("nz_DESI_coarse.txt",unpack=True); 
nzf=interp1d(zDESI_in,60*60*nzDESI_in,bounds_error=False,fill_value=0)
zDESI=zDESI_in[-1]*np.arange(256)/255.
nzDESI=nzf(zDESI)
bzDESI=0.84*pcs.growth_factor(1.)/np.array([pcs.growth_factor(1./(1+z)) for z in zDESI])
np.savetxt("nz_DESI.txt",np.transpose([zDESI,nzDESI]))
np.savetxt("bz_DESI.txt",np.transpose([zDESI,bzDESI]))
