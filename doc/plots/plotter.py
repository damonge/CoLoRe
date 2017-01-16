import numpy as np
import matplotlib.pyplot as plt

def mk_grid(nside,ax,col) :
    for i in np.arange(nside+1) :
        x0=(i+0.0)/nside;
        ax.plot([x0,x0],[0,1],col+'-')
        ax.plot([0,1],[x0,x0],col+'-')

def mk_grid_sphere(nr,ax,col) :
    for i in np.arange(nr) :
        rad=0.5*0.95*(i+1.)/nr
        ax.add_artist(plt.Circle((0.5,0.5),rad,color=col,fill=False,linewidth=2))
        nside=1; siz=2*np.pi*rad/nside; agrid=0.5/nr
        while siz>agrid :
            nside*=2; siz=2*np.pi*rad/nside; agrid=0.5/nr
        print i, agrid, siz, nside
        for j in np.arange(nside) :
            alpha=2*j*np.pi/nside
            rad0=0.5*0.95*(i+0.0)/nr
            radf=0.5*0.95*(i+1.0)/nr
            m=np.tan(alpha)
            x0=rad0*np.cos(alpha); xf=radf*np.cos(alpha)
            y0=m*x0; yf=m*xf;
            x0+=0.5; xf+=0.5; y0+=0.5; yf+=0.5
            ax.plot([x0,xf],[y0,yf],col+'-',linewidth=2)
    return nside
fig=plt.figure(figsize=[10,10])
ax=fig.add_subplot(111,aspect='equal')
mk_grid(16,ax,'k')
nside=mk_grid_sphere(8,ax,'r')
#nside=64
alpha=np.pi/2-2*4.5*np.pi/nside;
ax.arrow(0.5+0.5*0.95*np.cos(alpha),0.5+0.5*0.95*np.sin(alpha),
         -(0.5*0.95-0.04)*np.cos(alpha),-(0.5*0.95-0.04)*np.sin(alpha),
         head_width=0.015,head_length=0.04,fc='b',ec='b')
ax.tick_params(axis=u'both', which=u'both',length=0,labelsize=0)
plt.show()
