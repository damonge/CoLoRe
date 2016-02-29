import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad

ngals=1E7
zmax=0.9
#nz=100
#zarr=zmax*(np.arange(nz)+0.5)/nz
nz=101
zarr=zmax*np.arange(nz)/(nz-1.)
barr=1.*np.ones(nz)
narr=(zarr/0.2)**2*np.exp(-(zarr/0.2)**1.5)
nzi=interp1d(zarr,narr,bounds_error=False,fill_value=0)
nz2i=None

def integ_n(z,which) :
    if which==1 :
        nz=nzi(z)
    else :
        nz=nz2i(z)
    return 4*np.pi*(180/np.pi)**2*nz
norm=quad(integ_n,0,zmax,args=(1))[0]
narr*=ngals/norm
nz2i=interp1d(zarr,narr,bounds_error=False,fill_value=0)
print quad(integ_n,0,zmax,args=(2))[0]

plt.plot(zarr,narr); plt.show()

np.savetxt("nz_test.txt",np.transpose([zarr,narr]))
np.savetxt("bias_test.txt",np.transpose([zarr,barr]))
