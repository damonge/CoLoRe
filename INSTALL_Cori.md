# How to install CoLoRe in Cori (at NERSC)

## load relevant modules

You can choose whether to use GNU or Intel modules, here we use GNU

You can load the modules by hand every time, or setup a script with the 
following lines:

module swap PrgEnv-intel PrgEnv-gnu
module load gsl
module load fftw/3.3.4.11
module load cfitsio


## Compile libconfig

Clone the github repository for libconfig, and type:

mkdir $HOME/Install.cori
cd libconfig
autoreconf
./configure --prefix=$HOME/Install.cori
make 
make install


## Compile libsharp

Clone the github repository for libsharp, and type:

cd libsharp
autoreconf
./configure 
make 


