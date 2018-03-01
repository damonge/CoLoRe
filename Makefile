########## User-definable stuff ##########
#
###Compiler and compilation options
COMP_SER = cc
COMP_MPI = cc
OPTIONS = -Wall -O3 -std=c99
#
### Behavioural flags
#Use double precision integer (enable in general)
DEFINEFLAGS += -D_LONGIDS
#Use normalized bias model
#DEFINEFLAGS += -D_BIAS_MODEL_2
#Use linear bias model
#DEFINEFLAGS += -D_BIAS_MODEL_3
#Generate debug help. Only useful for development
DEFINEFLAGS += -D_DEBUG
#Use double precision floating point? Set to "yes" or "no"
USE_SINGLE_PRECISION = yes
#Add random perturbations to kappa from redshifts outside the box
ADD_EXTRA_KAPPA = yes
#Compile with HDF5 capability? Set to "yes" or "no"
USE_HDF5 = yes
#Use OMP parallelization? Set to "yes" or "no"
USE_OMP = yes
#Use MPI parallelization? Set to "yes" or "no"
USE_MPI = yes
#
###Path to libraries and headers
###If two or more of the dependencies reside in the same paths, only
###one instance is necessary.
#GSL
#GSL_INC = -I/add/path
#GSL_LIB = -L/add/path
GSL_INC = -I/usr/common/software/gsl/2.1/intel/include
GSL_LIB = -L/usr/common/software/gsl/2.1/intel/lib
#FFTW
FFTW_INC = -I/opt/cray/pe/fftw/3.3.6.2/haswell/include
FFTW_LIB = -L/opt/cray/pe/fftw/3.3.6.2/haswell/lib
#cfitsio
FITS_INC = -I/usr/common/software/cfitsio/3.370-reentrant/hsw/intel/include
FITS_LIB = -L/usr/common/software/cfitsio/3.370-reentrant/hsw/intel/lib
#cfitsio
HDF5_INC =
HDF5_LIB =
#libconfig
CONF_INC = -I/global/homes/d/damonge/include
CONF_LIB = -L/global/homes/d/damonge/lib
#healpix
HPIX_INC =
HPIX_LIB =
#libsharp
SHT_INC =
SHT_LIB =
#
########## End of user-definable ##########

USE_FITS = yes
DEFINEFLAGS += -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF 

ifeq ($(strip $(USE_OMP)),yes)
OPTIONS += -fopenmp
DEFINEFLAGS += -D_HAVE_OMP
endif #OMP

ifeq ($(strip $(USE_MPI)),yes)
DEFINEFLAGS += -D_HAVE_MPI
COMP_PAR = $(COMP_MPI)
else #MPI
COMP_PAR = $(COMP_SER)
endif #MPI

ifeq ($(strip $(USE_SINGLE_PRECISION)),yes)
DEFINEFLAGS += -D_SPREC

ifeq ($(strip $(USE_OMP)),yes)
LIB_FFTW += -lfftw3f_omp
endif #OMP
ifeq ($(strip $(USE_MPI)),yes)
LIB_FFTW += -lfftw3f_mpi
endif #MPI
LIB_FFTW += -lfftw3f

else #SINGLE_PRECISION

ifeq ($(strip $(USE_OMP)),yes)
LIB_FFTW += -lfftw3_omp
endif #OMP
ifeq ($(strip $(USE_MPI)),yes)
LIB_FFTW += -lfftw3_mpi 
endif #MPI

endif #SINGLE_PRECISION
LIB_FFTW += -lfftw3

# for fftlog
LIB_FFTW += -lfftw3

OPTIONS += $(DEFINEFLAGS)

INC_ALL = -I./src $(GSL_INC) $(FFTW_INC) $(FITS_INC) $(HDF5_INC) $(CONF_INC) $(SHT_INC) $(HPIX_INC)
LIB_ALL = $(GSL_LIB) $(FFTW_LIB) $(FITS_LIB) $(HDF5_LIB) $(CONF_LIB) $(SHT_LIB) $(HPIX_LIB) -lconfig -lgsl -lgslcblas $(LIB_FFTW) -lcfitsio -lchealpix
ifeq ($(strip $(ADD_EXTRA_KAPPA)),yes)
DEFINEFLAGS += -D_ADD_EXTRA_KAPPA -D_WITH_SHT
LIB_ALL += -lsharp -lfftpack -lc_utils
endif #EXTRA_KAPPA
ifeq ($(strip $(USE_HDF5)),yes)
DEFINEFLAGS += -D_HAVE_HDF5
LIB_ALL += -lhdf5 -lhdf5_hl -lz
endif #HDF5
ifeq ($(strip $(USE_FITS)),yes)
DEFINEFLAGS += -D_HAVE_FITS
#LIB_ALL += -lcfitsio
endif #FITS
LIB_ALL += -lm

COMMONO = src/common.o
COSMOMADO = src/cosmo_mad.o
COSMOO = src/cosmo.o
FOURIERO = src/fourier.o
DENSO = src/density.o
SRCSO = src/srcs.o
IMAPO = src/imap.o
KAPPAO = src/kappa.o
ISWO = src/isw.o
IOO = src/io.o
HPIXO = src/healpix_extra.o
BEAMO = src/beaming.o
PREDICTO = src/predictions.o
FFTLOGO = src/fftlog.o
MAIN = src/main.c
OFILES = $(COMMONO) $(COSMOMADO) $(COSMOO) $(FOURIERO) $(DENSO) $(BEAMO) $(IOO) $(HPIXO) $(SRCSO) $(IMAPO) $(KAPPAO) $(ISWO) $(FFTLOGO) $(PREDICTO)
#OFILES = $(COMMONO) $(COSMOMADO) $(COSMOO) $(FOURIERO) $(DENSO) $(LCO) $(IOO) $(HPIXO) $(PIXO) $(PREDICTO) $(FFTLOGO)

EXEC = CoLoRe

default : $(EXEC)

%.o : %.c
	$(COMP_CC) $(OPTIONS) $(INC_ALL) -c $< -o $@

$(OFILES) : COMP_CC := $(COMP_PAR)

$(EXEC) : $(OFILES)
	$(COMP_PAR) $(OPTIONS) $(INC_ALL) $(OFILES) $(MAIN) -o $(EXEC) $(LIB_ALL)

clean :
	rm -f src/*.o

cleaner :
	rm -f *~ src/*.o src/*~  $(EXEC)
