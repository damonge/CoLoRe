import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad

n1gals=1E6
n2gals=1E6
zmax=0.9
#nz=100
#zarr=zmax*(np.arange(nz)+0.5)/nz
nz=101
zarr=zmax*np.arange(nz)/(nz-1.)

b1arr=1.*np.ones(nz)
n1arr=(zarr/0.1)**2*np.exp(-(zarr/0.1)**1.5)
n1zi=interp1d(zarr,n1arr,bounds_error=False,fill_value=0)
n1z2i=None

b2arr=1.4*np.ones(nz)
n2arr=(zarr/0.25)**2*np.exp(-(zarr/0.25)**1.5)
n2zi=interp1d(zarr,n2arr,bounds_error=False,fill_value=0)
n2z2i=None

def integ_n1(z,which) :
    if which==1 :
        nz=n1zi(z)
    else :
        nz=n1z2i(z)
    return 4*np.pi*(180/np.pi)**2*nz

def integ_n2(z,which) :
    if which==1 :
        nz=n2zi(z)
    else :
        nz=n2z2i(z)
    return 4*np.pi*(180/np.pi)**2*nz

norm1=quad(integ_n1,0,zmax,args=(1))[0]
n1arr*=n1gals/norm1
n1z2i=interp1d(zarr,n1arr,bounds_error=False,fill_value=0)
norm2=quad(integ_n2,0,zmax,args=(1))[0]
n2arr*=n2gals/norm2
n2z2i=interp1d(zarr,n2arr,bounds_error=False,fill_value=0)

print quad(integ_n1,0,zmax,args=(2))[0], quad(integ_n2,0,zmax,args=(2))[0]

plt.plot(zarr,n1arr);
plt.plot(zarr,n2arr);
plt.show()

np.savetxt("nz1_test.txt",np.transpose([zarr,n1arr]))
np.savetxt("bias1_test.txt",np.transpose([zarr,b1arr]))
np.savetxt("nz2_test.txt",np.transpose([zarr,n2arr]))
np.savetxt("bias2_test.txt",np.transpose([zarr,b2arr]))
