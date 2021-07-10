import numpy as np
from scipy.integrate import simps


zs = np.linspace(0, 1., 1024)
z0 = 0.3
kz1 = np.exp(-5*(zs/z0))
kz2 = 4*kz1*(1-kz1)
bz1 = 1+zs
bz2 = np.sqrt(1+zs)
zx = 0.07
ng = 1E6
nz = (zs/zx)**2*np.exp(-(zs/zx)**1.5)
norm = simps(nz, x=zs)
nz *= ng/(norm*4*np.pi*(180/np.pi)**2)

np.savetxt("nz.txt", np.transpose([zs, nz]))
np.savetxt("kz1.txt", np.transpose([zs, kz1]))
np.savetxt("kz2.txt", np.transpose([zs, kz2]))
np.savetxt("bz1.txt", np.transpose([zs, bz1]))
np.savetxt("bz2.txt", np.transpose([zs, bz2]))