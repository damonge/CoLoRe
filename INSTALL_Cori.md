# How to install CoLoRe in Cori (at NERSC)

These instructions worked in November 1st 2019. 

## Load relevant modules

You can choose whether to use GNU or Intel modules, here we use GNU

You can load the modules by hand every time, or setup a script with the 
following lines:

```
module swap PrgEnv-intel PrgEnv-gnu
module load gsl
module load cray-fftw
module load cfitsio
module load openmpi
```

## Compile libconfig

Clone the github repository for libconfig: https://github.com/hyperrealm/libconfig, and type:
```
mkdir $HOME/Install.cori
cd libconfig
autoreconf
./configure --prefix=$HOME/Install.cori
make 
make install
```

## Compile libsharp

Clone the github repository for libsharp: https://github.com/Libsharp/libsharp, and type:
```
cd libsharp
autoreconf
./configure 
make 
```

## Modify the Makefile of CoLoRe 

Instructions below are to run Lyman alpha skewers in DESI.
```
COMP_SER = cc
COMP_MPI = cc
OPTIONS = -Wall -O3 -std=c99

DEFINEFLAGS += -D_LONGIDS
DEFINEFLAGS += -D_BIAS_MODEL_3 #
DEFINEFLAGS += -D_DEBUG
## Using single precision reduces the memory footprint
USE_SINGLE_PRECISION = yes 

ADD_EXTRA_KAPPA = no
USE_HDF5 = no
USE_FITS = yes
USE_OMP = yes
USE_MPI = yes

GSL_INC = -I${GSL_DIR}/include
GSL_LIB = -L${GSL_DIR}/lib

FFTW_INC = -I${FFTW_DIR}/include
FFTW_LIB = -L${FFTW_DIR}/lib

FITS_INC = -I${CFITSIO_DIR}/include
FITS_LIB = -L${CFITSIO_DIR}/lib

HPIX_INC = -I/global/common/cori/contrib/hpcosmo/hpcports_gnu-4.0/healpix-3.30.1_62c0405b-4.0/include
HPIX_LIB = -L/global/common/cori/contrib/hpcosmo/hpcports_gnu-4.0/healpix-3.30.1_62c0405b-4.0/lib

CONF_INC = -I/PATH/TO/LIBCONFIG/include
CONF_LIB = -L/PATH/TO/LIBCONFIG/lib

SHT_INC = -I/PATH/TO/LIBSHARP/auto/include
SHT_LIB = -L/PATH/TO/LIBSHARP/auto/lib
```

Make sure that the libraries used to compile `CoLoRe` are in `LD_LIBRARY_PATH` so they can be found at runtime:
`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${GSL_DIR}/lib:${FFTW_DIR}/lib:${CFITSIO_DIR}/lib:/PATH/TO/LIBCONFIG/lib:/PATH/TO/LIBSHARP/auto/lib`, etc.


## Possible issue on compute nodes

Although `CoLoRe` should run well on the login nodes after following these instructions, it can sometimes crash on the compute nodes. All instances of these have been solved by calling `module unload craype-hugepages2M` before running the code.
