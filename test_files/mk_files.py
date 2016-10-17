#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

nz=256
z0=0
zf=1.0
zx=0.10
pow1=1.5
ngal=1.4E6
zx2=0.60
ngal2=2e6
pow2=4.0

nux=1420.
n_nu=15
z0nu=0.1
zfnu=0.5
nu0=nux/(1+zfnu)
nuf=nux/(1+z0nu)
dnu=(nuf-nu0)/n_nu
nu0_arr=nu0+dnu*np.arange(n_nu)
nuf_arr=nu0+dnu*(1+np.arange(n_nu))

zarr=z0+(zf-z0)*np.arange(nz)/(nz-1.)
nzarr=(zarr/zx)**2*np.exp(-(zarr/zx)**pow1); 
nzarr2=(zarr/zx2)**2*np.exp(-(zarr/zx2)**pow2); 
bzarr=np.ones_like(zarr)
bzarr2=2*np.ones_like(zarr)
tzarr=5.5919E-2+2.3242E-1*zarr-2.4136E-2*zarr**2.
norm=ngal/(4*np.pi*(180/np.pi)**2*np.sum(nzarr)*(zf-z0)/nz)
nzarr*=norm
norm2=ngal2/(4*np.pi*(180/np.pi)**2*np.sum(nzarr2)*(zf-z0)/nz)
nzarr2*=norm2

np.savetxt("Nz_test.txt",np.transpose([zarr,nzarr]))
np.savetxt("Nz2_test.txt",np.transpose([zarr,nzarr2]))
np.savetxt("Tz_test.txt",np.transpose([zarr,tzarr]))
np.savetxt("Bz_test.txt",np.transpose([zarr,bzarr]))
np.savetxt("Bz2_test.txt",np.transpose([zarr,bzarr2]))
np.savetxt("nuTable.txt",np.transpose([nu0_arr,nuf_arr]))

plt.plot(zarr,nzarr); plt.plot(zarr,nzarr2); plt.show()
plt.plot(zarr,bzarr); plt.show()
plt.plot(zarr,tzarr); plt.show()
