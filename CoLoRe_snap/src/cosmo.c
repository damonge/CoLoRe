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
#include "common.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

static void int_error_handle(int status,double result,
                             double error)
{
  //////
  // Error handler for gsl
  if(isnan(result)) {
    fprintf(stderr,"CoLoRe: NAN found \n");
    exit(1);
  }
  else{
    if(status==GSL_EROUND)
      fprintf(stderr,"CoLoRe: Roundoff error: %lE %lE \n",result,error);
    else if(status==GSL_EMAXITER)
      fprintf(stderr,"CoLoRe: Ran out of iterations: %lE %lE \n",result,error);
    else if(status==GSL_ESING)
      fprintf(stderr,"CoLoRe: Singularity found: %lE %lE \n",result,error);
    else if(status==GSL_EDIVERGE) {
      fprintf(stderr,"CoLoRe: Integral seems to diverge: %lE %lE \n",
	      result,error);
    }
    else if(status==GSL_ETOL) {
      fprintf(stderr,"CoLoRe: Can't reach tolerance: %lE %lE : %lE %%\n",
              result,error,100*error/result);
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

static double wind(double x,int setwf)
{
  //////
  // Window function:
  //  setwf=0 -> top hat
  //  setwf=1 -> gaussian
  //  setwf=2 -> sharp-k
  if(setwf==1) { //Gaussian
    return exp(-0.5*x*x);
  }
  else if(setwf==0) { //TopHat
    if(x<0.1) {
      return 1.-0.1*x*x+0.003571429*x*x*x*x-6.61376E-5*x*x*x*x*x*x
        +7.51563E-7*x*x*x*x*x*x*x*x;
    }
    else
      return 3*(sin(x)-x*cos(x))/(x*x*x);
  }
  else if(setwf==2) { //Sharp-k
    if(x<1) return 1;
    else return 0;
  }
  else
    return -1;
}

double pk_linear0(ParamCoLoRe *par,double lgk)
{
  //////
  // Linear power spectrum at redshift 0.
  // Extrapolated to ~k^ns for small k and
  // to k^{-3} for large k
  double pk;
  int ik=(int)((lgk-par->logkmin)*par->idlogk);
  
  if(ik<0)
    pk=par->pkarr[0]*pow(10,par->n_scal*(lgk-par->logkmin));
  else if(ik<par->numk)
    pk=par->pkarr[ik]+(lgk-par->logkarr[ik])*(par->pkarr[ik+1]-par->pkarr[ik])*par->idlogk;
  else
    pk=par->pkarr[par->numk-1]*pow(10,-3*(lgk-par->logkmax));

  return pk;
}

static double j_bessel_0(double x)
{
  //////
  // Bessel's j0 function
  if(x>0.1)
    return sin(x)/x;
  else
    return 1.-0.166667*x*x+0.00833333*x*x*x*x-
      0.000198413*x*x*x*x*x*x+2.75573E-6*x*x*x*x*x*x*x*x;
}

typedef struct { //Param struct for integrals
  double r;
  double R1;
  double R2;
  int wf1;
  int wf2;
  ParamCoLoRe *par;
} xiparam;

static double integxiL_O(double kk,void *params)
{
  //////
  // integrand for xi_L (for oscillatory integration)
  double dum;
  double x1,x2,xr;
  xiparam *par;
  double lgk=log10(kk);
  par=(xiparam *)params;

  x1=kk*(par->R1);
  x2=kk*(par->R2);
  xr=kk*(par->r);

  dum=TWOPIPIINV*pk_linear0(par->par,lgk)*kk*kk*
    wind(x1,par->wf1)*wind(x2,par->wf2)/xr;

  return dum;
}

static double integxiL_NO(double logk,void *params)
{
  //////
  // integrand for xi_L (including oscillatory term in j0)
  double dum;
  double x1,x2,xr;
  xiparam *par;
  par=(xiparam *)params;

  double kk=pow(10,logk);
  x1=kk*(par->R1);
  x2=kk*(par->R2);
  xr=kk*(par->r);

  dum=TWOPIPIINVLOGTEN*pk_linear0(par->par,logk)*kk*kk*kk*
    wind(x1,par->wf1)*wind(x2,par->wf2)*j_bessel_0(xr);

  return dum;
}

static double xi2p_L(ParamCoLoRe *par,double r,double R1,double R2,
		     char *wf1,char *wf2,double errfac)
{
  //////
  // Correlation function between the linear density contrast smoothed
  // with window function (wf1,R1) and with window function (wf2,R2)
  // at two points separated by a distance r:
  //              <delta_(R1,wf1)(x)*delta_(R2,wf2)(x+r)>
  gsl_function integrand;
  double relerrt=1E-4;
  double integral,errintegral;
  xiparam xpar;
  double lim=MIN(R1,R2);
  lim/=r;

  xpar.r=r;
  xpar.R1=R1;
  xpar.R2=R2;
  xpar.par=par;
  if(!strcmp(wf1,"Gauss"))
    xpar.wf1=1;
  else if(!strcmp(wf1,"TopHat"))
    xpar.wf1=0;
  else {
    fprintf(stderr,"CoLoRe: Unknown window function %s \n",wf1);
    exit(1);
  }
  if(!strcmp(wf2,"Gauss"))
    xpar.wf2=1;
  else if(!strcmp(wf2,"TopHat"))
    xpar.wf2=0;
  else {
    fprintf(stderr,"CoLoRe: Unknown window function %s \n",wf2);
    exit(1);
  }

  gsl_integration_workspace *w
    =gsl_integration_workspace_alloc(1000);
  integrand.params=&xpar;
  if(lim>=1) {
    integrand.function=&integxiL_NO;
    int stat=gsl_integration_qagil(&integrand,par->logkmax,0,relerrt,1000,w,
                                   &integral,&errintegral);
    int_error_handle(stat,integral,errintegral);
  }
  else {
    lim*=errfac;
    gsl_integration_workspace *cw
      =gsl_integration_workspace_alloc(1000);
    gsl_integration_qawo_table *wf
      =gsl_integration_qawo_table_alloc(r,0.1,GSL_INTEG_SINE,100);

    integrand.function=&integxiL_O;
    int stat=gsl_integration_qawf(&integrand,0,relerrt*lim,1000,
                           w,cw,wf,&integral,&errintegral);
    int_error_handle(stat,integral,errintegral);

    gsl_integration_qawo_table_free(wf);
    gsl_integration_workspace_free(cw);
  }
  gsl_integration_workspace_free(w);

  return integral;
}

static double sigL2(ParamCoLoRe *par,double R1,double R2,char *wf1,char *wf2)
{
  //////
  // Covariance between the linear density contrast smoothed with
  // window function (wf1,R1) and with window function (wf2,R2) at
  // the same point:  <delta_(R1,wf1)(x)*delta_(R2,wf2)(x)>
  return xi2p_L(par,0,R1,R2,wf1,wf2,1);
}

static void pk_linear_set(ParamCoLoRe *par)
{
  //////
  // Reads linear power spectrum. CAMB format expected.
  int ii;
  double kk,ppk;
  FILE *fpk;

  print_info("Reading P_k from file: %s\n",par->fnamePk);
  fpk=fopen(par->fnamePk,"r");
  if(fpk==NULL) error_open_file(par->fnamePk);
  par->numk=linecount(fpk);
  par->logkarr=(double *)my_malloc(par->numk*sizeof(double));
  par->pkarr=(double *)my_malloc(par->numk*sizeof(double));
  rewind(fpk);
  for(ii=0;ii<par->numk;ii++) {
    int stat=fscanf(fpk,"%lf %lf",&kk,&ppk);
    if(stat!=2) error_read_line(par->fnamePk,ii+1);
    par->pkarr[ii]=ppk;
    par->logkarr[ii]=log10(kk); //log(k) in h Mpc^-1
  }
  fclose(fpk);
  
  par->logkmin=par->logkarr[0];
  par->logkmax=par->logkarr[par->numk-1];
  par->idlogk=(par->numk-1)/(par->logkmax-par->logkmin);

  //Re-interpolate just in case the file is not equi-spaced in log10(k)
  gsl_interp_accel *intacc=gsl_interp_accel_alloc();
  gsl_spline *spline=gsl_spline_alloc(gsl_interp_cspline,par->numk);
  gsl_spline_init(spline,par->logkarr,par->pkarr,par->numk);
  for(ii=0;ii<par->numk-1;ii++) {
    double lk=par->logkmin+ii/(par->idlogk);
    par->pkarr[ii]=gsl_spline_eval(spline,lk,intacc);
    par->logkarr[ii]=lk;
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(intacc);

  // normalize
  double norm_pk=par->sig8*par->sig8/sigL2(par,8,8,"TopHat","TopHat");
  print_info("  Original sigma8=%lf\n",
	     sqrt(sigL2(par,8,8,"TopHat","TopHat")));
  for(ii=0;ii<par->numk;ii++)
    par->pkarr[ii]*=norm_pk;

  double r_effective=sqrt(par->r2_smooth+pow(0.45*par->l_box/par->n_grid,2));
  par->sigma2_gauss=sigL2(par,r_effective,r_effective,"Gauss","Gauss");
#ifdef _DEBUG
  print_info("  Sigma_Gauss should be %lf\n",sqrt(par->sigma2_gauss));
#endif //_DEBUG
}

void cosmo_set(ParamCoLoRe *par)
{
  //////
  // This initializes the cosmological model
  // at redshift z_s
  Csm_params *pars=csm_params_new();
  csm_unset_gsl_eh();
  csm_background_set(pars,par->OmegaM,par->OmegaL,par->OmegaB,par->weos,0,par->hhub,2.275);
  double h0=csm_hubble(pars, 1);
  double fgrowth0=csm_f_growth(pars, 1);
  double growth0=csm_growth_factor(pars,1);
  par->prefac_lensing = 1.5*h0*h0*par->OmegaM;

  // Compute growth factor at this redshift
  double a=1/(1+par->z_snap);
  double om=csm_omega_m(pars,a);
  double d1=csm_growth_factor(pars,a)/growth0;
  double d2=-0.42857142857*d1*d1*pow(om,-0.00699300699);
  double fz=csm_f_growth(pars,a);
  double hz=csm_hubble(pars,a);
  par->growth_d1=d1;
  par->growth_d2=d2;
  par->growth_dv=(d1*hz*fz)/(fgrowth0*h0); //This is for the comoving velocity
  par->ihub=1/hz;
  par->fgrowth_0=fgrowth0;
  par->hubble_0=h0;

  pk_linear_set(par);

  csm_params_free(pars);
  //FREEE SLPINES!
}
