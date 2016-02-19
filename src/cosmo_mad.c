///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CoLoRe.                                      //
//                                                                   //
// CoLoRe is free software: you can redistribute it and/or modify    //
// it under the terms of the GNU General Public License as published //
// by the Free Software Foundation, either version 3 of the License, //
// or (at your option) any later version.                            //
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
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "common.h"
#include "cosmo_mad.h"

//Number of points for the a(t) relation
#define CSM_NINTERP_A 1000

/**** Verbosity ****/
static int csm_flag_verbose=1;

/**** Error handler ****/
static gsl_error_handler_t *csm_gsl_error_handler_old;


/****************************/
/*     General routines     */
/****************************/
static void int_error_handle(int status,double result,
                             double error)
{
  //////
  // Error handler for gsl
  if(isnan(result)) {
    fprintf(stderr,"CoLoRe: NAN found \n");
  }
  else{
    if(status==GSL_EROUND) {
      fprintf(stderr,"CoLoRe: Roundoff error: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_EMAXITER) {
      fprintf(stderr,"CoLoRe: Ran out of iterations: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_ESING) {
      fprintf(stderr,"CoLoRe: Singularity found: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_EDIVERGE) {
      fprintf(stderr,"CoLoRe: Integral seems to diverge: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_ETOL) {
      fprintf(stderr,"CoLoRe: Can't reach tolerance: %lE %lE\n",
	      result,error);
    }
    else if(status==GSL_EUNDRFLW)
      fprintf(stderr,"CoLoRe: Underflow: %lE %lE \n",result,error);
    else if(status==GSL_EDOM) {
      fprintf(stderr,"CoLoRe: Outside interpolation range!! %lE %lE\n",
	      result,error);
      exit(1);
    }
    else if(status) {
      fprintf(stderr,"CoLoRe: Unknown error code %d %lf %lf \n",
	      status,result,error);
      exit(1);
    }
  }
}

void csm_unset_gsl_eh(void)
{
  //////
  // Disables GSL default error handler
  csm_gsl_error_handler_old=gsl_set_error_handler_off();
}

void csm_set_verbosity(int verb)
{
  //////
  // Sets verbosity level
  csm_flag_verbose=verb;
}

static void csm_bg_params_free(Csm_bg_params *par)
{
  //////
  // bg_params destructor
  free(par);
}

static Csm_bg_params *csm_bg_params_new(void)
{
  //////
  // bg_params creator
  // Default Planck parameters
  Csm_bg_params *bg=(Csm_bg_params *)malloc(sizeof(Csm_bg_params));
  if(bg==NULL) {
    fprintf(stderr,"CoLoRe: Out of memory!\n");
    exit(1);
  }
  bg->OM=0.315;
  bg->OL=0.685;
  bg->OB=0.049;
  bg->h=0.673;
  bg->w0=-1;
  bg->wa=0;
  bg->TCMB=2.725;
  bg->OK=0;
  bg->ksign=0;
  bg->normalDE=1;
  bg->constantw=1;

  return bg;
}

void csm_params_free(Csm_params *par)
{
  //////
  // csm_params destructor
  if(par->bg_params_set)
    csm_bg_params_free(par->bg);
  par->bg_params_set=0;
  free(par);
}

Csm_params *csm_params_new(void)
{
  //////
  // csm_params destructor
  Csm_params *par=(Csm_params *)malloc(sizeof(Csm_params));
  if(par==NULL) {
    fprintf(stderr,"CoLoRe: Out of memory!\n");
    exit(1);
  }
  par->bg_params_set=0;
  par->bg=NULL;

  return par;
}

/****************************/
/*   Background cosmology   */
/****************************/

static double aeqmin(Csm_bg_params *par)
{
  //////
  // Returns MIN(aeq_k,aeq_L), where aeq_i is the scale
  // factor at equality of M with k or L.
  double aeqK=1;
  double aeqL=1;
  
  if(par->ksign!=0)
    aeqK=par->OM/fabs(par->OK);

  if (par->OL!=0) {
    if(par->normalDE)
      aeqL=pow(par->OM/par->OL,0.333);
    else
      aeqL=pow(par->OM/par->OL,-1/(3*par->w0));
  }

  return CSM_MIN(aeqK,aeqL);
}

static double sinn(double x,int sign)
{
  //////
  //         { sin(x)  , if k==1
  // sinn(x)={  x      , if k==0
  //         { sinh(x) , if k==-1
  double dum;

  if(sign==-1)
    dum=sinh(x);
  else if(sign==1)
    dum=sin(x);
  else
    dum=x;
  
  return dum;
}

static double naHm1(double a,void *params)
{
  //////
  // H0/(a*H[a])
  double dum;
  Csm_bg_params *par=(Csm_bg_params *)params;

  if(par->normalDE) 
    dum=sqrt(a/(par->OM+par->OL*a*a*a+par->OK*a));
  else if(par->constantw) {
    dum=sqrt(a/(par->OM+par->OL*pow(a,-3*par->w0)+par->OK*a));
  }
  else {
    dum=sqrt(a/(par->OM+par->OL*pow(a,-3*(par->w0+par->wa))*
		exp(3*(a-1)*par->wa)+par->OK*a));
  }

  return dum;
}

static double na2Hm1(double a,void *params)
{
  //////
  // H0/(a^2*H[a])
  double dum;
  Csm_bg_params *par=(Csm_bg_params *)params;

  if(par->normalDE) 
    dum=1/sqrt(a*(par->OM+par->OL*a*a*a+par->OK*a));
  else if(par->constantw) {
    dum=1/sqrt(a*(par->OM+par->OL*pow(a,-3*par->w0)+par->OK*a));
  }
  else {
    dum=1/sqrt(a*(par->OM+par->OL*pow(a,-3*(par->w0+par->wa))*
		  exp(3*(a-1)*par->wa)+par->OK*a));
  }

  return dum;
}

static double naHm3(double a,void *params)
{
  //////
  // (H0/(a*H[a]))^3
  double dum;
  Csm_bg_params *par=(Csm_bg_params *)params;

  if(par->normalDE) {
    dum=sqrt(a/(par->OM+par->OL*a*a*a+par->OK*a));
  }
  else if(par->constantw) {
    dum=sqrt(a/(par->OM+par->OL*pow(a,-3*par->w0)+par->OK*a));
  }
  else {
    dum=sqrt(a/(par->OM+par->OL*pow(a,-3*(par->w0+par->wa))*
		  exp(3*(a-1)*par->wa)+par->OK*a));
  }

  return dum*dum*dum;
}

static void parthor(Csm_bg_params *par,
		    double aa,double *ph,double *delph)
{
  //////
  // Particle horizon. The returned value is in Mpc/h.
  double alim;

  alim=0.01*par->a_equality;
  if(aa<=alim) {
    *ph=2*sqrt(aa/par->OM)*CSM_HMPC;
    *delph=0;
  }
  else {
    double relerrt=1E-6;
    double integral,errintegral,int0;
    size_t sdum;
    gsl_function integrand;
    int stat;

    int0=2*sqrt(alim/par->OM);
    integrand.function=&na2Hm1;
    integrand.params=par;

    stat=gsl_integration_qng(&integrand,alim,aa,0,relerrt,
			     &integral,&errintegral,&sdum);
    int_error_handle(stat,integral,errintegral);
    *ph=(int0+integral)*CSM_HMPC;
    //    *delph=errintegral*CSM_HMPC;
  }
}

static void gfac(Csm_bg_params *par,
		 double aa,double *gf,double *delgf)
{
  //////
  // Growth factor, normalized to gfac(a<<a_eq)~a
  double alim,int0;

  alim=0.01*par->a_equality;
  int0=0.4*sqrt(alim*alim*alim*alim*alim/(par->OM*par->OM*par->OM));
  if(aa<=alim) {
    *gf=aa;
    *delgf=0;
  }
  else {
    double relerrt=1E-4;
    double integral,errintegral;
    size_t sdum;
    gsl_function integrand;
    int stat;
   
    integrand.function=&naHm3;
    integrand.params=par;
    
    stat=gsl_integration_qng(&integrand,alim,aa,0,relerrt,
			     &integral,&errintegral,&sdum);
    int_error_handle(stat,integral,errintegral);
    *gf=(int0+integral)*2.5*par->OM/(aa*naHm1(aa,par));
  }
}

double csm_hubble(Csm_params *par,double aa)
{
  //////
  // Hubble rate at aa in h/Mpc
  if(par->bg->normalDE) {
    return sqrt((par->bg->OM+par->bg->OL*aa*aa*aa+par->bg->OK*aa)/
		(aa*aa*aa))/CSM_HMPC;
  }
  else if(par->bg->constantw) {
    return sqrt((par->bg->OM+par->bg->OL*pow(aa,-3*par->bg->w0)+
		 par->bg->OK*aa)/(aa*aa*aa))/CSM_HMPC;
  }
  else {
    return sqrt((par->bg->OM+par->bg->OL*pow(aa,-3*(par->bg->w0+par->bg->wa))*
		 exp(3*par->bg->wa*(aa-1))+par->bg->OK*aa)/(aa*aa*aa))/CSM_HMPC;
  }
}

double csm_particle_horizon(Csm_params *par,double aa)
{
  //////
  // Particle horizon
  double ph, eph;

  parthor(par->bg,aa,&ph,&eph);
  return ph;
}

double csm_radial_comoving_distance(Csm_params *par,double aa)
{
  //////
  // chi(a)
  double rcd;

  rcd=csm_particle_horizon(par,aa);
  return par->bg->phorizon-rcd;
}

double csm_curvature_comoving_distance(Csm_params *par,double aa)
{
  //////
  // r(a)
  if(par->bg->ksign==0)
    return csm_radial_comoving_distance(par,aa);
  else {
    double dum;
    double ksq=sqrt(fabs(par->bg->OK));
    dum=csm_radial_comoving_distance(par,aa)/CSM_HMPC;
    dum=sinn(ksq*dum,par->bg->ksign)/ksq;
    return dum*CSM_HMPC;
  }
}

double csm_growth_factor(Csm_params *par,double aa)
{
  //////
  // D(a)
  double gf,egf;
  
  gfac(par->bg,aa,&gf,&egf);

  return gf;
}

double csm_f_growth(Csm_params *par,double aa)
{
  //////
  // f(a)= d ln(D)/d ln(a)
  // Using the integral definition for D(a)
  double apow,coeff,Da;

  Da=csm_growth_factor(par,aa);
  if(par->bg->normalDE) {
    coeff=0;
    apow=aa*aa*aa;
  }
  else if(par->bg->constantw) {
    coeff=(1+par->bg->w0);
    apow=pow(aa,-3*par->bg->w0);
  }
  else {
    coeff=(1+par->bg->w0+par->bg->wa*(1-aa));
    apow=pow(aa,-3*(par->bg->w0+par->bg->wa))*exp(3*par->bg->wa*(aa-1));
  }

  return 0.5*(5*par->bg->OM*aa/Da-
	      (3*par->bg->OM+3*coeff*par->bg->OL*apow+2*par->bg->OK*aa))/
    (par->bg->OM+par->bg->OL*apow+par->bg->OK*aa);
}

void csm_background_set(Csm_params *par,
			double OmegaM,double OmegaL,double OmegaB,
			double ww,double wwa,double hh,double T_CMB)
{
  //////
  // This initializes the cosmological parameters.
  if(par->bg_params_set)
    csm_bg_params_free(par->bg);
  par->bg_params_set=0;

  par->bg=csm_bg_params_new();
  par->bg_params_set=1;
  par->bg->h=hh;
  par->bg->w0=ww;
  par->bg->wa=wwa;
  par->bg->OM=OmegaM;
  par->bg->OL=OmegaL;
  par->bg->OK=1-par->bg->OM-par->bg->OL;
  par->bg->OB=OmegaB;
  par->bg->TCMB=T_CMB;

  //Check parameters
  if(fabs(par->bg->wa)<1E-6) {
    if(fabs(par->bg->w0+1)<1E-6) {
      par->bg->constantw=1;
      par->bg->normalDE=1;
    }
    else {
      par->bg->constantw=1;
      par->bg->normalDE=0;
    }
  }
  else {
    par->bg->constantw=0;
    par->bg->normalDE=0;
  }

  if(fabs(par->bg->OK)<1E-6)
    par->bg->ksign=0;
  else if(par->bg->OK>0)
    par->bg->ksign=-1;
  else
    par->bg->ksign=1; 

  if(par->bg->OM<=0) {
    fprintf(stderr,"CoLoRe: Wrong matter parameter %.3lf \n",
	    par->bg->OM);
    exit(1);
  }
  if(par->bg->OM<par->bg->OB) {
    fprintf(stderr,"CoLoRe: Wrong M/B parameter %.3lf > %.3lf \n",
	    par->bg->OB,par->bg->OM);
    exit(1);
  }
  if(par->bg->w0>-0.333333) {
    fprintf(stderr,"CoLoRe: DE is too exotic (w=%.3lf \n",par->bg->w0);
    exit(1);
  }
  if(par->bg->TCMB<0) {
    fprintf(stderr,"CoLoRe: Wrong CMB temperature %.3lf \n",
	    par->bg->TCMB);
    exit(1);
  }

  if(csm_flag_verbose) {
    print_info("The cosmological model is:\n");
    print_info(" O_M=%.3f O_L=%.3f O_K=%.3f\n",
	       par->bg->OM,par->bg->OL,par->bg->OK);
    print_info(" O_B=%.3f w=%.3f h=%.3f\n",
	       par->bg->OB,par->bg->w0,par->bg->h);
    if(par->bg->ksign==0)
      print_info(" Flat universe, ");
    else if(par->bg->ksign==1)
      print_info(" Closed universe, ");
    else if(par->bg->ksign==-1)
      print_info(" Open universe, ");
    if(par->bg->normalDE)
      print_info("standard cosmological constant\n");
    else {
      print_info("non-standard dark energy");
      if(par->bg->constantw)
	print_info("\n");
      else
	print_info(": w(a) = %.3lf + %.3lf*(1-a) \n",par->bg->w0,par->bg->wa);
    }
  }
  par->bg->a_equality=aeqmin(par->bg);
  par->bg->phorizon=csm_particle_horizon(par,1);
  par->bg->growth0=csm_growth_factor(par,1);
  if(csm_flag_verbose) {
    print_info("\n Time of equality: a_eq=%.5lf\n",par->bg->a_equality);
    print_info(" Particle horizon: ");
    print_info("chi_H(0)=%.3lE Mpc/h\n",par->bg->phorizon);
    print_info(" Present growth factor: ");
    print_info("D_0=%.3lf\n\n",par->bg->growth0);
  }
}
