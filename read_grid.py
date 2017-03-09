import numpy as np
import sys as sys
import matplotlib.pyplot as plt
from matplotlib import cm

def read_grid(prefix) :
    f=open(prefix+"_0.dat","rb")
    nfiles,size_float=np.fromfile(f,dtype=np.int32,count=2)
    lbox=np.fromfile(f,dtype=np.float64,count=1)[0]
    ngrid=np.fromfile(f,dtype=np.int32,count=1)[0]
    f.close()

    if size_float==4 :
        f_type=np.float32
    else :
        f_type=np.float64

    grid_out=np.zeros([3,ngrid,ngrid])
    for ifil in np.arange(nfiles) :
        f=open(prefix+"_%d.dat"%ifil,"rb")
        nf,sz=np.fromfile(f,dtype=np.int32,count=2)
        lb=np.fromfile(f,dtype=np.float64,count=1)
        ng,nz_here,iz0_here=np.fromfile(f,dtype=np.int32,count=3)
        for iz in np.arange(nz_here) :
            print ifil,iz
            g=np.fromfile(f,dtype=f_type,count=ng*ng).reshape([ng,ng])
            if iz0_here+iz==ngrid/2 :
                grid_out[2,:,:]=g
            grid_out[0,:,iz0_here+iz]=g[ngrid/2,:]
            grid_out[1,:,iz0_here+iz]=g[:,ngrid/2]
        f.close()

    return ngrid,lbox,np.array(grid_out)

if len(sys.argv)!= 2 :
    print "Usage: read_grid.py prefix"
    exit(1)
prefix=sys.argv[1]

ng,lb,dens=read_grid(prefix)
print "Ngrid=%d"%ng
print "Lbox=%.3lf Mpc/h"%lb
def plot_slice(slic) :
    plt.figure()
    plt.imshow(slic,origin='lower',interpolation='nearest');
    plt.colorbar()
    plt.savefig("sl%d.png"%ib,bbox_inches='tigh')
for ib in [0,1,2] :
    plot_slice(dens[ib,:,:],ib) 
#plot_slice(dens[1,:,:])
#plot_slice(dens[2,:,:])
#plt.show()
