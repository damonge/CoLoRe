# How to install CoLoRe in Cori (at NERSC)

These instructions worked in November 2nd 2017. 

## Load relevant modules

You can choose whether to use GNU or Intel modules, here we use GNU

You can load the modules by hand every time, or setup a script with the 
following lines:

```
module swap PrgEnv-intel PrgEnv-gnu
module load gsl
module load fftw/3.3.4.11
module load cfitsio
```

## Compile libconfig

Clone the github repository for libconfig, and type:
```
mkdir $HOME/Install.cori
cd libconfig
autoreconf
./configure --prefix=$HOME/Install.cori
make 
make install
```

## Compile libsharp

Clone the github repository for libsharp, and type:
```
cd libsharp
autoreconf
./configure 
make 
```

## Modify the Makefile of CoLoRe 

Instructions below are to run Lyman alpha skewers in DESI.
```
COMP_SER = gcc
COMP_MPI = mpicc
OPTIONS = -Wall -O3 -std=c99

DEFINEFLAGS += -D_LONGIDS
DEFINEFLAGS += -D_BIAS_MODEL_3
DEFINEFLAGS += -D_DEBUG

USE_SINGLE_PRECISION = yes

ADD_EXTRA_KAPPA = no
USE_HDF5 = no
USE_FITS = yes
USE_OMP = yes
USE_MPI = yes

GSL_INC = -I/usr/common/software/gsl/2.1/gnu/include
GSL_LIB = -L/usr/common/software/gsl/2.1/gnu/lib -lgsl -lgslcblas

FFTW_INC = -I/opt/cray/pe/fftw/3.3.4.11/haswell/include/
FFTW_LIB = -L/opt/cray/pe/fftw/3.3.4.11/haswell/lib/

FITS_INC = -I/usr/common/software/cfitsio/3.370-reentrant/hsw/gnu/include
FITS_LIB = -L/usr/common/software/cfitsio/3.370-reentrant/hsw/gnu/lib -lcfitsio

HPIX_INC = -I/global/common/cori/contrib/hpcosmo/hpcports_gnu-4.0/healpix-3.30.1_62c0405b-4.0/include
HPIX_LIB = -L//global/common/cori/contrib/hpcosmo/hpcports_gnu-4.0/healpix-3.30.1_62c0405b-4.0/lib

CONF_INC = -I/global/homes/f/font/Install.cori/include
CONF_LIB = -L/global/homes/f/font/Install.cori/lib

SHT_INC = -I/global/homes/f/font/Programs/Others/libsharp/auto/include
SHT_LIB = -L/global/homes/f/font/Programs/Others/libsharp/auto/lib
```
