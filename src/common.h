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

//Interpolation type
#ifndef INTERP_TYPE_SKW
#define INTERP_TYPE_SKW INTERP_CIC
#endif //INTERP_TYPE_SKW
#ifndef INTERP_TYPE_SHEAR
#define INTERP_TYPE_SHEAR INTERP_NGP
#endif //INTERP_TYPE_SHEAR

//dr_par = dx/NSAMP_RAD
#ifndef NSAMP_RAD
#define NSAMP_RAD 1
#endif //NSAMP_RAD

//Maximum allowed healpix resolution
#define NSIDE_MAX_HPX 8192

//Background tags
#define BG_Z 1000
#define BG_D1 1001
#define BG_D2 1002
#define BG_V1 1003
#define BG_PD 1004
#define BG_IH 1005
#define BG_NZ_SRCS 1006
#define BG_BZ_SRCS 1007
#define BG_NORM_SRCS 1008
#define BG_TZ_IMAP 1009
#define BG_BZ_IMAP 1010
#define BG_NORM_IMAP 1011

//Density field parameters
#define DZ_SIGMA 0.05
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
#define NA 5001
#define NPOP_MAX 10
#define NPLANES_MAX 100

#ifdef _HAVE_MPI
#ifdef _SPREC
#define FLOUBLE_MPI MPI_FLOAT
#else //_SPREC
#define FLOUBLE_MPI MPI_DOUBLE
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
  int *ipix;
} CatalogCartesian;

typedef struct {
  float ra;     //Right ascension
  float dec;    //Declination
  float z0;     //Cosmological redshift
  float dz_rsd; //RSD contribution
  float e1;
  float e2;
} Src;

typedef struct {
  int nsrc;
  Src *srcs;
  int nr;
  double rmax;
  double dr;
  double idr;
  int has_skw;
  float *d_skw;
  float *v_skw;
} Catalog;

typedef struct {
  int nside; //Resolution parameter
  long num_pix;
  long *listpix;
  double *pos;
  int nr; //Number of spherical shells
  flouble *r0; //r_min in this shell
  flouble *rf; //r_max in this shell
  flouble *data;
  int *nadd;
} HealpixShells;

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
  double fgrowth_0; //Growth rate at z=0
  double hubble_0; //Expansion rate at z=0 (inverse length units)
  double prefac_lensing; //3*O_M*H_0^2/2
  double z_max; //Maximum redshift
  double z_min; //Minimum redshift
  double r_max; //Maximum radial comoving distance
  double r_min; //Minimum radial comoving distance
  double a_arr_a2r[NA]; //Array of redshifts used to compute r(z)
  double r_arr_a2r[NA]; //Array of comoving distances used to compute r(z)
  double z_arr_r2z[NA]; //Array of redshifts used to compute z(r)
  double r_arr_r2z[NA]; //Array of comoving distances used to compute z(r), D_d(r), D_v(r), 1/H(z)
  double growth_d_arr[NA]; //Array of density growth factors used to compute D_d(r)
  double growth_d2_arr[NA]; //Array of density growth factors used to compute D_d(r)
  double growth_v_arr[NA]; //Array of velocity growth factors used to compute D_v(r)
  double growth_pd_arr[NA]; //Array of potential derivative factors used to compute \dot{\phi}
  double ihub_arr[NA]; //Array of 1/H(z)
  double glob_idr; //1/dr, where dr is the radial comoving distance interval used in the arrays above
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
  flouble l_box; //Box size for the cartesian grid
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
  flouble *slice_left; //Dummy array to store grid cells coming from the left node
  flouble *slice_right; //Dummy array to store grid cells coming from the right node
  double sigma2_gauss; //Variance of the cartesian density field

  //IO parameters
  char prefixOut[256]; //Output prefix
  int output_format; //0-> ASCII, 1-> FITS, 2-> HDF5
  double pos_obs[3]; //Observer position

  //Tracers
  // Sources
  int do_srcs; //Do we include sources?
  int do_skewers; //Do we include skewer information?
  int do_lensing; //Do we need to compute the lensing potential?
  int n_srcs; //Number of source types
  char fnameBzSrcs[NPOP_MAX][256]; //Files containing b(z) for each source type
  char fnameNzSrcs[NPOP_MAX][256]; //Files containing dN/dzdOmega (in deg^-2)
  double *srcs_nz_arr[NPOP_MAX];
  double *srcs_bz_arr[NPOP_MAX];
  double *srcs_norm_arr[NPOP_MAX];
  double norm_srcs_0[NPOP_MAX]; //Bottom edge of spline for density normalization
  double norm_srcs_f[NPOP_MAX]; //Top edge of spline for density normalization  
  int shear_srcs[NPOP_MAX]; //Do we do lensing for this source type?
  int skw_srcs[NPOP_MAX]; //Do we want to store skewers for each source of this type?
  long *nsources_c_this; //Number of sources initially found in this node
  CatalogCartesian **cats_c; //Galaxy positions initially stored in this node
  long *nsources_this; //Number of sources finally found in this node
  Catalog **cats; //Final galaxy properties
  // Intensity mapping
  int do_imap; //Do we include intensity mapping
  int n_imap; //Number of IM species
  char fnameBzImap[NPOP_MAX][256]; //Files containing b(z) for each IM species
  char fnameTzImap[NPOP_MAX][256]; //Files containing T(z) for each IM species
  char fnameNuImap[NPOP_MAX][256]; //Files containing frequency table for each IM species
  double *imap_tz_arr[NPOP_MAX];
  double *imap_bz_arr[NPOP_MAX];
  double *imap_norm_arr[NPOP_MAX];
  double norm_imap_0[NPOP_MAX]; //Bottom edge of spline for density normalization
  double norm_imap_f[NPOP_MAX]; //Top edge of spline for density normalization
  int nside_imap[NPOP_MAX]; //Output angular resolution for each IM species
  double nu0_imap[NPOP_MAX]; //Rest-frame frequency for each IM species
  HealpixShells **imap; //intensity maps for each IM species
  // Kappa
  int do_kappa; //Do you want to output kappa maps?
  int n_kappa; //How many maps?
  double z_kappa_out[NPLANES_MAX]; //Array of source plane redshifts
  int nside_kappa;
  HealpixShells *kmap; //Kappa maps at each redshift
#ifdef _ADD_EXTRA_KAPPA
  int *need_extra_kappa;
  flouble **fl_mean_extra_kappa;
  flouble **cl_extra_kappa;
#endif //_ADD_EXTRA_KAPPA
  // ISW
  int do_isw; //Do you want to output isw maps?
  int n_isw; //How many maps?
  double z_isw_out[NPLANES_MAX]; //Array of source plane redshifts
  int nside_isw;
  HealpixShells *pd_map; //Isw maps at each redshift
#ifdef _ADD_EXTRA_KAPPA
  int *need_extra_isw;
  flouble **fl_mean_extra_isw;
  flouble **cl_extra_isw;
#endif //_ADD_EXTRA_KAPPA

  //Lognormal predictions
  int do_pred;
  int just_do_pred;
  double pred_dz;

  //Beam distribution
  int nside_base; //Minimum n_side used in the pixelization
  int npix_base; //Corresponding number of pixels
  int need_beaming; //Do we need spherical voxels at all?
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
void get_radial_params(double rmax,int ngrid,int *nr,double *dr);
//void get_random_angles(gsl_rng *rng,int ipix_nest,int ipix0,int nside,double *th,double *phi);
HealpixShells *hp_shell_alloc(int nside,int nside_base,int nr);
void hp_shell_free(HealpixShells *shell);
CatalogCartesian *catalog_cartesian_alloc(int nsrcs);
void catalog_cartesian_free(CatalogCartesian *cat);
Catalog *catalog_alloc(int nsrcs,int has_skw,double rmax,int ng);
void catalog_free(Catalog *cat);

static inline double bias_model(double d,double b)
{
  if(d<=-1)
    return 0;
#ifdef _BIAS_MODEL_2
  return pow(1+d,b)/pow(1+d*d,0.5*(b-1));
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
double r_of_z(ParamCoLoRe *par,double z);
double get_bg(ParamCoLoRe *par,double r,int tag,int ipop);
void compute_tracer_cosmo(ParamCoLoRe *par);

//////
// Functions defined in io.c
ParamCoLoRe *read_run_params(char *fname,int test_memory);
void write_density_grid(ParamCoLoRe *par,char *prefix_dens);
void write_lpt(ParamCoLoRe *par,unsigned long long npart,flouble *x,flouble *y,flouble *z);
void write_srcs(ParamCoLoRe *par);
void write_imap(ParamCoLoRe *par);
void write_kappa(ParamCoLoRe *par);
void write_isw(ParamCoLoRe *par);
void param_colore_free(ParamCoLoRe *par);


/////
// Functions defined in predictions.c
void write_predictions(ParamCoLoRe *par);


//////
// Functions defined in fourier.c
void init_fftw(ParamCoLoRe *par);
void allocate_fftw(ParamCoLoRe *par);
void create_cartesian_fields(ParamCoLoRe *par);
void end_fftw(ParamCoLoRe *par);
void fftw_wrap_c2r(int ng,dftw_complex *pin,flouble *pout);
void fftw_wrap_r2c(int ng,flouble *pin,dftw_complex *pout);


//////
// Functions defined in beaming.c
int interpolate_from_grid(ParamCoLoRe *par,double *x,
			  flouble *d,flouble v[3],flouble t[6],flouble *pd,
			  int flag_return,int interp_type);
void get_beam_properties(ParamCoLoRe *par);


//////
// Functions defined in density.c
void compute_physical_density_field(ParamCoLoRe *par);
void compute_density_normalization(ParamCoLoRe *par);


//////
// Functions defined in srcs.c
void srcs_set_cartesian(ParamCoLoRe *par);
void srcs_distribute(ParamCoLoRe *par);
void srcs_get_local_properties(ParamCoLoRe *par);
void srcs_beams_preproc(ParamCoLoRe *par);
void srcs_get_beam_properties(ParamCoLoRe *par);
void srcs_beams_postproc(ParamCoLoRe *par);


//////
// Functions defined in imap.c
void imap_set_cartesian(ParamCoLoRe *par);
void imap_distribute(ParamCoLoRe *par);
void imap_get_local_properties(ParamCoLoRe *par);
void imap_beams_preproc(ParamCoLoRe *par);
void imap_get_beam_properties(ParamCoLoRe *par);
void imap_beams_postproc(ParamCoLoRe *par);


//////
// Functions defined in kappa.c
void kappa_set_cartesian(ParamCoLoRe *par);
void kappa_distribute(ParamCoLoRe *par);
void kappa_get_local_properties(ParamCoLoRe *par);
void kappa_beams_preproc(ParamCoLoRe *par);
void kappa_get_beam_properties(ParamCoLoRe *par);
void kappa_beams_postproc(ParamCoLoRe *par);


//////
// Functions defined in isw.c
void isw_set_cartesian(ParamCoLoRe *par);
void isw_distribute(ParamCoLoRe *par);
void isw_get_local_properties(ParamCoLoRe *par);
void isw_beams_preproc(ParamCoLoRe *par);
void isw_get_beam_properties(ParamCoLoRe *par);
void isw_beams_postproc(ParamCoLoRe *par);


//////
// Defined in healpix_extra.c
long he_nside2npix(long nside);
double he_pixel_area(int nside);
long he_ang2pix(long nside,double cth,double phi);
void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname,int isnest);
flouble *he_read_healpix_map(char *fname,long *nside,int nfield);
int he_ring_num(long nside,double z);
long *he_query_strip(long nside,double theta1,double theta2,long *npix_strip);
void he_ring2nest_inplace(flouble *map_in,long nside);
void he_nest2ring_inplace(flouble *map_in,long nside);
void he_udgrade(flouble *map_in,long nside_in,flouble *map_out,long nside_out,int nest);
#ifdef _WITH_SHT
#ifdef _SPREC
#define SHT_TYPE 0
#else //_SPREC
#define SHT_TYPE SHARP_DP
#endif //_SPREC
#define HE_MAX_SHT 32
#define HE_FWHM2SIGMA 0.00012352884853326381 //Transforms FWHM in arcmin to sigma_G in rad:
long he_nalms(int lmax);
long he_indexlm(int l,int m,int lmax);
void he_alm2map(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
void he_map2alm(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms);
void he_alm2cl(fcomplex **alms_1,fcomplex **alms_2,
	       int nmaps_1,int nmaps_2,int pol_1,int pol_2,flouble **cls,int lmax);
void he_anafast(flouble **maps_1,flouble **maps_2,
		int nmaps_1,int nmaps_2,int pol_1,int pol_2,
		flouble **cls,int nside,int lmax);
flouble *he_generate_beam_window(int lmax,flouble fwhm_amin);
void he_alter_alm(int lmax,flouble fwhm_amin,fcomplex *alm_in,
		  fcomplex *alm_out,flouble *window);
flouble *he_synfast(flouble *cl,int nside,int lmax,unsigned int seed);
//HE_NT
#ifdef _WITH_NEEDLET
#define HE_NBAND_NX 512
#define HE_NORM_FT 2.2522836206907617
#define HE_NL_INTPREC 1E-6
#define HE_NT_NSIDE_MIN 32
typedef struct {
  double b;
  double inv_b;
  gsl_spline *b_spline;
  gsl_interp_accel *b_intacc;
  int nside0;
  int nj;
  int *nside_arr;
  int *lmax_arr;
  flouble **b_arr;
} HE_nt_param;
void he_nt_end(HE_nt_param *par);
HE_nt_param *he_nt_init(flouble b_nt,int nside0);
flouble ***he_alloc_needlet(HE_nt_param *par,int pol);
void he_free_needlet(HE_nt_param *par,int pol,flouble ***nt);
void he_nt_get_window(HE_nt_param *par,int j,flouble *b);
fcomplex **he_map2needlet(HE_nt_param *par,flouble **map,flouble ***nt,
			  int return_alm,int pol,int qu_in,int qu_out);
fcomplex **he_needlet2map(HE_nt_param *par,flouble **map,flouble ***nt,
			  int return_alm,int pol,int qu_in,int qu_out);
#endif //_WITH_NEEDLET
#endif //_WITH_SHT


#endif //_COMMON_
