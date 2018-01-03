import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

hdulist = fits.open('out_srcs_s0_0.fits')

#First HDU contains the source catalog
print hdulist[1].header.keys
plt.figure(); plt.hist(hdulist[1].data['Z_COSMO'],bins=100)
print " "

#Second HDU contains the density skewers as a FITS image
#The skewers have the same ordering as the sources in the catalog
#(i.e. skewer hdulist[2].data[i,:] corresponds to source hdulist[1].data[i])
id=np.argmax(hdulist[1].data['Z_COSMO'])
print hdulist[2].header.keys
plt.figure(); plt.plot(hdulist[4].data['R'],hdulist[2].data[id]);
plt.xlabel('$r\\,\\,[{\\rm Mpc}/h]$',fontsize=18);
plt.ylabel('$\\delta$',fontsize=18)
print " "

#Third HDU contains the velocity skewers. The units of the velocity are
#such that the skewers contain the redshift distortion associated with 
#the peculiar velocity field
print hdulist[3].header.keys
plt.figure(); plt.plot(hdulist[4].data['R'],hdulist[3].data[id]);
plt.xlabel('$r\\,\\,[{\\rm Mpc}/h]$',fontsize=18);
plt.ylabel('$\\delta z_{\\rm RSD}$',fontsize=18)
print " "

#Fourth HDU is a table containing background cosmological quantities at
#the distances where the skewers are sampled (see the use of
#hdulist[4].data['R'] in the previous examples
print hdulist[4].header.keys
plt.figure();
plt.plot(hdulist[4].data['Z'],hdulist[4].data['R']*0.001,label='$r(z)\\,[{\\rm Gpc}/h]$')
plt.plot(hdulist[4].data['Z'],hdulist[4].data['D'],label='$D_\\delta(z)$')
plt.plot(hdulist[4].data['Z'],hdulist[4].data['V'],label='$D_v(z)$')
plt.legend(loc='lower right')
plt.xlabel('$z$',fontsize=18)
print " "

plt.show()
