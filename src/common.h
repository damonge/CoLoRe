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
#define PT_CEA 0
#define PT_CAR 1
#define PT_HPX 2

#define INTERP_NGP 0
#define INTERP_CIC 1

//Pixelization type
#ifndef PIXTYPE
#define PIXTYPE PT_HPX
#endif //PIXTYPE

//Interpolation type
#ifndef INTERP_TYPE
#define INTERP_TYPE INTERP_CIC
#endif //INTERP_TYPE

//Resolution parameter for nearest onion shell
#ifndef NSIDE_ONION_BASE
#define NSIDE_ONION_BASE 2
#endif //NSIDE_ONION_BASE

//dr_par = FAC_CART2SPH_PAR * dx
#ifndef FAC_CART2SPH_PAR
#define FAC_CART2SPH_PAR 1.
#endif //FAC_CART2SPH_PAR

//dr_perp = FAC_CART2SPH_PERP * dx
#ifndef FAC_CART2SPH_PERP
#define FAC_CART2SPH_PERP 1.
#endif //FAC_CART2SPH_PERP

//#sub-voxel divisions in r
#ifndef NSUB_PAR
#define NSUB_PAR 1
#endif //NSUB_PAR

//#sub-voxel divisions in r
#ifndef NSUB_PERP
#define NSUB_PERP 1
#endif //NSUB_PERP

//sqrt(#random points per pixel for IM)
#ifndef NSUB_IMAP_PERP
#define NSUB_IMAP_PERP 4
#endif //NSUB_IMAP_PERP

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

#ifdef _HAVE_MPI

#ifdef _SPREC
#define FLOUBLE_MPI MPI_FLOAT
#else //_SPREC
#define FLOUBLE_MPI MPI_DOUBLE
#endif //_SPREC

#ifdef _LONGIDS
#define LINT_MPI MPI_INT
#else //_LONGIDS
#define LINT_MPI MPI_LONG
#endif //_LONGIDS

#endif //_HAVE_MPI

#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
typedef int lint;
#endif //_LONGIDS

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

typedef struct {
  float ra;     //Right ascension
  float dec;    //Declination
  float z0;     //Cosmological redshift
  float dz_rsd; //RSD contribution
  float e1;
  float e2;
} Src;

typedef struct {
  int nr;
  flouble *r0_arr;
  flouble *rf_arr;
  int *nside_arr;
  int *nside_ratio_arr;
  int *iphi0_arr;
  int *icth0_arr;
  int *num_pix;
} OnionInfo;

typedef struct {
  int nside; //Resolution parameter
  long num_pix;
  long *listpix;
  int nr; //Number of spherical shells
  flouble *r0; //r_min in this shell
  flouble *rf; //r_max in this shell
  flouble *data;
  int *nadd;
} HealpixShells;

typedef struct {
  char fnamePk[256]; //File containing power spectrum
  double OmegaM; //Cosmological parameters
  double OmegaL; //Cosmological parameters
  double OmegaB; //Cosmological parameters
  double hhub; //Cosmological parameters
  double weos; //Cosmological parameters
  double n_scal; //Cosmological parameters
  double sig8; //Cosmological parameters
  //Derived parameters
  double fgrowth_0; //Growth rate at z=0
  double hubble_0; //Expansion rate at z=0 (inverse length units)
  double prefac_lensing; //3*O_M*H_0^2/2
  double z_max; //Maximum redshift
  double z_min; //Minimum redshift
  double r_max; //Maximum radial comoving distance
  double r_min; //Minimum radial comoving distance
  double r2_smooth; //Square of the smoothing scale
  int smooth_potential; //Do we smooth the newtonian potential as well?
  int do_smoothing; //Are we smoothing the density field?

#ifdef _DEBUG
  FILE *f_dbg; //File into which all debug info is written
#endif //_DEBUG

  //Only used in common.c
  int numk; //Number of k-values
  double logkmax; //Maximum log10(k)
  double logkmin; //Minimum log10(k)
  double idlogk; //1/D(log10(k))
  double *logkarr; //Array of log10(k) values (units of h/Mpc)
  double *pkarr; //Array of power spectrum values (units of (Mpc/h)^3)
  double a_arr_a2r[NA]; //Array of redshifts used to compute r(z)
  double r_arr_a2r[NA]; //Array of comoving distances used to compute r(z)
  double z_arr_r2z[NA]; //Array of redshifts used to compute z(r)
  double r_arr_r2z[NA]; //Array of comoving distances used to compute z(r), D_d(r), D_v(r), 1/H(z)
  double growth_d_arr[NA]; //Array of density growth factors used to compute D_d(r)
  double growth_v_arr[NA]; //Array of velocity growth factors used to compute D_v(r)
  double growth_pd_arr[NA]; //Array of potential derivative factors used to compute \dot{\phi}
  double ihub_arr[NA]; //Array of 1/H(z)
  double glob_idr; //1/dr, where dr is the radial comoving distance interval used in the arrays above

  unsigned int seed_rng; //RNG seed

  int n_grid; //Number of cells per side for the Cartesian grid
  double l_box; //Box size for the cartesian grid
  int nz_here; //Number of cells in the z-direction stored in this node
  int iz0_here; //index of the first cell in the z-direction stored in this node
  int nz_max;
  int *nz_all;
  int *iz0_all;

  char prefixOut[256]; //Output prefix
  int output_format; //0-> ASCII, 1-> FITS, 2-> HDF5
  int output_density; //Do you want to output the density grid?
  double pos_obs[3]; //Observer position

  dftw_complex *grid_dens_f; //Fourier-space grid for the density field
  flouble *grid_dens; //Real-space grid for the density field
  dftw_complex *grid_npot_f; //Fourier-space grid for the Newtonian potential
  flouble *grid_npot; //Real-space grid for the Newtonian potential
  flouble *slice_left; //Dummy array to store grid cells coming from the left node
  flouble *slice_right; //Dummy array to store grid cells coming from the right node

  double sigma2_gauss; //Variance of the cartesian density field
  double z0_sigma2; //Start z for sigma(z)
  double zf_sigma2; //End z for sigma(z)
  double sigma2_0; //Start sigma(z)
  double sigma2_f; //End sigma(z)
  gsl_interp_accel *intacc_sigma2_z; //Splines sigma(z)
  gsl_spline *spline_sigma2_z; //Splines sigma(z)

  int need_onions; //Do we need spherical voxels at all?
  int do_lensing; //Do we need to compute the lensing potential?
  int nside_base; //Minimum n_side used in the pixelization
  int n_beams_here; //Number of beams stored in this node for the lightcone
  OnionInfo **oi_beams; //Onion beams stored in this node
  flouble ***dens_beams; //Density beams
  flouble ***vrad_beams; //v_r beams
  flouble ***p_xx_beams; //phi_xx beams
  flouble ***p_xy_beams; //phi_xy beams
  flouble ***p_yy_beams; //phi_yy beams
  flouble ***pdot_beams; //phi_t beams
  int ***nsrc_beams; //Beams with total number of sources

  int do_sources; //Do we include sources
  int n_srcs; //Number of source types
  char fnameBzSrcs[NPOP_MAX][256]; //Files containing b(z) for each source type
  char fnameNzSrcs[NPOP_MAX][256]; //Files containing dN/dzdOmega 
                                   //(number of sources per redshift interval and square degree)
  gsl_spline *spline_srcs_bz[NPOP_MAX]; //Spline for b(z)
  gsl_spline *spline_srcs_nz[NPOP_MAX]; //Spline for n(z)
  gsl_interp_accel *intacc_srcs[NPOP_MAX]; //Spline accelerator for sources
  int shear_srcs[NPOP_MAX]; //Do we do lensing for this source type?
  lint *nsources_this; //Number of sources found in this node
  Src **srcs; //Galaxy objects stored in this node

  int do_imap; //Do we include intensity mapping
  int n_imap; //Number of IM species
  char fnameBzImap[NPOP_MAX][256]; //Files containing b(z) for each IM species
  char fnameTzImap[NPOP_MAX][256]; //Files containing T(z) for each IM species
  char fnameNuImap[NPOP_MAX][256]; //Files containing frequency table for each IM species
  gsl_spline *spline_imap_bz[NPOP_MAX]; //Spline for b(z)
  gsl_spline *spline_imap_tz[NPOP_MAX]; //Spline for T(z)
  gsl_interp_accel *intacc_imap[NPOP_MAX]; //Spline accelerator for IM
  int nside_imap[NPOP_MAX]; //Output angular resolution for each IM species
  double nu0_imap[NPOP_MAX]; //Rest-frame frequency for each IM species
  HealpixShells **imap; //intensity maps for each IM species

  int do_kappa; //Do you want to output kappa maps?
  int n_kappa; //How many maps?
  double z_kappa_out[NPOP_MAX]; //Array of source plane redshifts
  int nside_kappa;
  HealpixShells *kmap; //Kappa maps at each redshift
#ifdef _ADD_EXTRA_KAPPA
  int *need_extra_kappa;
  flouble **fl_mean_extra_kappa;
  flouble **cl_extra_kappa;
#endif //_ADD_EXTRA_KAPPA

  int do_isw; //Do you want to output isw maps?
  int n_isw; //How many maps?
  double z_isw_out[NPOP_MAX]; //Array of source plane redshifts
  int nside_isw;
  HealpixShells *pd_map; //Isw maps at each redshift
#ifdef _ADD_EXTRA_ISW
  int *need_extra_isw;
  flouble **fl_mean_extra_isw;
  flouble **cl_extra_isw;
#endif //_ADD_EXTRA_ISW

  int do_pred;
  double pred_dz;
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
void timer(int i);
gsl_rng *init_rng(unsigned int seed);
double rng_01(gsl_rng *rng);
int rng_poisson(double lambda,gsl_rng *rng);
void rng_delta_gauss(double *module,double *phase,
		     gsl_rng *rng,double sigma2);
void rng_gauss(gsl_rng *rng,double *r1,double *r2);
void end_rng(gsl_rng *rng);
OnionInfo **alloc_onion_info_beams(ParamCoLoRe *par);
void free_onion_info(OnionInfo *oi);
void alloc_beams(ParamCoLoRe *par);
void free_beams(ParamCoLoRe *par);
flouble get_res(int nside);
long get_npix(int nside);
void get_vec(int ipix_nest,int iphi_0,int icth_0,int nside,int nside_ratio,double *u);
void get_random_angles(gsl_rng *rng,int ipix_nest,int iphi_0,int icth_0,int nside,int nside_ratio,
		       double *th,double *phi);


//////
// Functions defined in cosmo.c
double pk_linear0(ParamCoLoRe *par,double lgk);
void cosmo_set(ParamCoLoRe *par);
double sigma2_of_z(ParamCoLoRe *par,double z);
double r_of_z(ParamCoLoRe *par,double z);
double z_of_r(ParamCoLoRe *par,double r);
double dgrowth_of_r(ParamCoLoRe *par,double r);
double vgrowth_of_r(ParamCoLoRe *par,double r);
double pdgrowth_of_r(ParamCoLoRe *par,double r);
double ihub_of_r(ParamCoLoRe *par,double r);
double ndens_of_z_srcs(ParamCoLoRe *par,double z,int ipop);
double bias_of_z_srcs(ParamCoLoRe *par,double z,int ipop);
double temp_of_z_imap(ParamCoLoRe *par,double z,int ipop);
double bias_of_z_imap(ParamCoLoRe *par,double z,int ipop);
void free_hp_shell(HealpixShells *shell);
HealpixShells *new_hp_shell(int nside,int nr);


//////
// Functions defined in io.c
ParamCoLoRe *read_run_params(char *fname);
void write_catalog(ParamCoLoRe *par);
void write_imap(ParamCoLoRe *par);
void write_kappa(ParamCoLoRe *par);
void write_isw(ParamCoLoRe *par);
void write_grids(ParamCoLoRe *par);
void param_colore_free(ParamCoLoRe *par);


/////
// Functions defined in predictions.h
void write_predictions(ParamCoLoRe *par);


//////
// Functions defined in fourier.c
void init_fftw(ParamCoLoRe *par);
void create_cartesian_fields(ParamCoLoRe *par);
void end_fftw(ParamCoLoRe *par);


//////
// Functions defined in pixelization.c
void pixelize(ParamCoLoRe *par);


//////
// Functions defined in grid_tools.c
void integrate_lensing(ParamCoLoRe *par);
void integrate_isw(ParamCoLoRe *par);
void get_sources(ParamCoLoRe *par);
void get_imap(ParamCoLoRe *par);
void get_kappa(ParamCoLoRe *par);
void get_isw(ParamCoLoRe *par);


//////
// Defined in healpix_extra.c
long he_nside2npix(long nside);
long he_ang2pix(long nside,double cth,double phi);
void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname);
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
