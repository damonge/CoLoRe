///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CoLoRe.                                      //
//                                                                   //
// CoLoRe is free software: you can redistribute it and/or modify it //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// CoLoRe is distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CoLoRe.  If not, see <http://www.gnu.org/licenses/>.   //
//                                                                   //
///////////////////////////////////////////////////////////////////////
#ifndef _COMMON_
#define _COMMON_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#ifdef _HAVE_OMP
#include <omp.h>
#endif //_HAVE_OMP
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "fftw3.h"
#ifdef _HAVE_MPI
#include <mpi.h>
#include <fftw3-mpi.h>
#endif //_HAVE_MPI
#include <libconfig.h>
#ifdef _HAVE_FITS
#include <fitsio.h>
#endif //_HAVE_FITS
#ifdef _HAVE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif //_HAVE_HDF5
#include <chealpix.h>
#ifdef _WITH_SHT
#include <sharp_almhelpers.h>
#include <sharp_geomhelpers.h>
#include <sharp.h>
#ifdef _WITH_NEEDLET
#include <gsl/gsl_integration.h>
#endif //_WITH_NEEDLET
#endif //_WITH_SHT
#include "cosmo_mad.h"

/////////
// Interpolation parameters

#define INTERP_NGP 0
#define INTERP_CIC 1
#define INTERP_TSC 2
#define RETURN_DENS 1
#define RETURN_VEL  2
#define RETURN_TID  4
#define RETURN_PDOT 8
#define RETURN_GAUSS 16

//Interpolation type
#ifndef INTERP_TYPE_SKW
#define INTERP_TYPE_SKW INTERP_CIC
#endif //INTERP_TYPE_SKW
#ifndef INTERP_TYPE_LENSING
#define INTERP_TYPE_LENSING INTERP_NGP
#endif //INTERP_TYPE_LENSING

//Density field parameters
#define DENS_TYPE_LGNR 0
#define DENS_TYPE_1LPT 1
#define DENS_TYPE_2LPT 2
#define DENS_TYPE_CLIP 3

// End of interpolation parameters
/////////

#define DYNAMIC_SIZE 1
#define RTOD 57.2957795
#define DTOR 0.01745329251
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

#define TWOPIPIINVLOGTEN  0.1166503235296796 //ln(10)/(2*pi^2)
#define TWOPIPIINV  0.05066059182116889 //1/(2*pi^2)
#define NPOP_MAX 10

#ifdef _HAVE_MPI
#ifdef _SPREC
#define FLOUBLE_MPI MPI_FLOAT
#define FCOMPLEX_MPI MPI_C_FLOAT_COMPLEX
#else //_SPREC
#define FLOUBLE_MPI MPI_DOUBLE
#define FCOMPLEX_MPI MPI_C_DOUBLE_COMPLEX
#endif //_SPREC
#endif //_HAVE_MPI

#ifdef _SPREC
typedef float flouble;
typedef float complex fcomplex;
typedef fftwf_complex dftw_complex;
#else //_SPREC
typedef double flouble;
typedef double complex fcomplex;
typedef fftw_complex dftw_complex;
#endif //_SPREC

//Defined in common.c
extern int NodeThis;
extern int NodeLeft;
extern int NodeRight;
extern int NNodes;
extern int IThread0;
extern int MPIThreadsOK;

#define NPOS_CC 4
typedef struct {
  int nsrc;
  float *pos;
} CatalogCartesian;

typedef struct {

#ifdef _DEBUG
  FILE *f_dbg; //File into which all debug info is written
#endif //_DEBUG

  //Cosmological parameters
  // Background
  double OmegaM; //Cosmological parameters
  double OmegaL; //Cosmological parameters
  double OmegaB; //Cosmological parameters
  double hhub; //Cosmological parameters
  double weos; //Cosmological parameters
  double n_scal; //Cosmological parameters
  double sig8; //Cosmological parameters
  double prefac_lensing; //3*O_M*H_0^2/2
  // Power spectra
  char fnamePk[256]; //File containing power spectrum
  int numk; //Number of k-values
  double logkmax; //Maximum log10(k)
  double logkmin; //Minimum log10(k)
  double idlogk; //1/D(log10(k))
  double *logkarr; //Array of log10(k) values (units of h/Mpc)
  double *pkarr; //Array of power spectrum values (units of (Mpc/h)^3)

  //Density parameters
  // Density methods
  int output_density; //Do you want to output the density grid?
  double r2_smooth; //Square of the smoothing scale
  int do_smoothing; //Are we smoothing the density field?
  int smooth_potential; //Do we smooth the newtonian potential as well?
  int dens_type; //Method to produce the density field
  int lpt_interp_type;
  double lpt_buffer_fraction; //Fraction of memory saved for buffer particles
  int output_lpt;
  unsigned int seed_rng; //RNG seed
  // Box parameters
  int n_grid; //Number of cells per side for the Cartesian grid
  double l_box; //Box size for the cartesian grid
  double z_snap; //Snapshot redshift
  flouble growth_d1; //Growth factor at z_snap
  flouble growth_d2; //2nd-orther growth at z_snap
  flouble growth_dv; //Velocity growth at z_snap
  flouble growth_fz; //Growth rate at z_snap
  flouble ihub; //Hubble scale (in Mpc/h)
  flouble fgrowth_0; //f(z=0)
  flouble hubble_0; //inverse hubble scale at z=0
  int nz_here; //Number of cells in the z-direction stored in this node
  int iz0_here; //index of the first cell in the z-direction stored in this node
  int nz_max;
  int *nz_all;
  int *iz0_all;
  double z0_norm;
  double zf_norm;
  // Density grids
  dftw_complex *grid_dens_f; //Fourier-space grid for the density field
  flouble *grid_dens; //Real-space grid for the density field
  dftw_complex *grid_npot_f; //Fourier-space grid for the Newtonian potential
  flouble *grid_npot; //Real-space grid for the Newtonian potential
  dftw_complex *grid_eta_f; //Fourier-space grid for the RSD field
  flouble *grid_eta; //Real-space grid for the RSD field
  flouble *slice_left; //Dummy array to store grid cells coming from the left node
  flouble *slice_right; //Dummy array to store grid cells coming from the right node
  dftw_complex *edge_planes; //Buffer used to populate k=0 and k=ngrid/2 planse
  double sigma2_gauss; //Variance of the cartesian density field

  //IO parameters
  char prefixOut[256]; //Output prefix
  int output_format; //0-> ASCII, 1-> FITS, 2-> HDF5

  //Tracers
  // Sources
  int do_srcs; //Do we include sources?
  int n_srcs; //Number of source types
  long *nsources_this; //Number of sources initially found in this node
  CatalogCartesian **cats; //Galaxy positions initially stored in this node
  double *bias;
  double *threshold;
  double *ndens;
  double *dens_norm;
} ParamCoLoRe;

void mpi_init(int* p_argc,char*** p_argv);
void *my_malloc(size_t size);
void *my_calloc(size_t nmemb,size_t size);
size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream);
void error_open_file(char *fname);
void error_read_line(char *fname,int nlin);
void print_info(char *fmt,...);
void report_error(int level,char *fmt,...);
int linecount(FILE *f);
int *ind_sort(int n,flouble *arr);
void timer(int i);
gsl_rng *init_rng(unsigned int seed);
double rng_01(gsl_rng *rng);
int rng_poisson(double lambda,gsl_rng *rng);
void rng_delta_gauss(double *module,double *phase,
		     gsl_rng *rng,double sigma2);
void rng_gauss(gsl_rng *rng,double *r1,double *r2);
void end_rng(gsl_rng *rng);
unsigned long long get_max_memory(ParamCoLoRe *par,int just_test);
CatalogCartesian *catalog_cartesian_alloc(int nsrcs);
void catalog_cartesian_free(CatalogCartesian *cat);

static inline double bias_model(double d,double b, double t)
{
  if(d<=-1 || d < t)
    return 0;
#ifdef _BIAS_MODEL_2
  if(d < 0)
    return exp(b*d/(1+d));
  else
    return 1+b*d;
#elif defined _BIAS_MODEL_3
  if(1+b*d>0)
    return 1+b*d;
  else
    return 0;
#else //_BIAS_MODEL
  return pow(1+d,b);
#endif //_BIAS_MODEL
}


//////
// Functions defined in cosmo.c
double pk_linear0(ParamCoLoRe *par,double lgk);
void cosmo_set(ParamCoLoRe *par);

//////
// Functions defined in io.c
ParamCoLoRe *read_run_params(char *fname,int test_memory);
void write_density_grid(ParamCoLoRe *par,char *prefix_dens);
void write_eta_grid(ParamCoLoRe *par);
void write_lpt(ParamCoLoRe *par,unsigned long long npart,flouble *x,flouble *y,flouble *z);
void write_srcs(ParamCoLoRe *par);
void param_colore_free(ParamCoLoRe *par);

//////
// Functions defined in fourier.c
void init_fftw(ParamCoLoRe *par);
void allocate_fftw(ParamCoLoRe *par);
void create_cartesian_fields(ParamCoLoRe *par);
void end_fftw(ParamCoLoRe *par);
void fftw_wrap_c2r(int ng,dftw_complex *pin,flouble *pout);
void fftw_wrap_r2c(int ng,flouble *pin,dftw_complex *pout);

//////
// Functions defined in density.c
void compute_physical_density_field(ParamCoLoRe *par);
void compute_density_normalization(ParamCoLoRe *par);

//////
// Functions defined in srcs.c
void srcs_set_cartesian(ParamCoLoRe *par);

#endif //_COMMON_
