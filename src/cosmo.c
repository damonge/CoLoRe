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

#define CL_TYPE_KAPPA 0
#define CL_TYPE_ISW 1

static double f_of_r_linear(ParamCoLoRe *par,double r,double *f,double f0,double ff)
{
  if(r<=0) return f0;
  else if(r>=par->r_arr_r2z[NA-1]) return ff;
  else {
    int ir=(int)(r*par->glob_idr);
    return  f[ir]+(f[ir+1]-f[ir])*(r-par->r_arr_r2z[ir])*par->glob_idr;
  }
}

double get_bg(ParamCoLoRe *par,double r,int tag,int ipop)
{
  if(tag==BG_Z)
    return f_of_r_linear(par,r,par->z_arr_r2z,0,par->z_arr_r2z[NA-1]);
  else if(tag==BG_D1)
    return f_of_r_linear(par,r,par->growth_d_arr,1,par->growth_d_arr[NA-1]);
  else if(tag==BG_D2) {
    return f_of_r_linear(par,r,par->growth_d2_arr,par->growth_d2_arr[0],
			 par->growth_d2_arr[NA-1]);
  }
  else if(tag==BG_V1)
    return f_of_r_linear(par,r,par->growth_v_arr,par->growth_v_arr[0],par->growth_v_arr[NA-1]);
  else if(tag==BG_PD) {
    return f_of_r_linear(par,r,par->growth_pd_arr,par->growth_pd_arr[0],
			 par->growth_pd_arr[NA-1]);
  }
  else if(tag==BG_IH)
    return f_of_r_linear(par,r,par->ihub_arr,par->ihub_arr[0],par->ihub_arr[NA-1]);
  else if(tag==BG_NZ_SRCS)
    return f_of_r_linear(par,r,par->srcs_nz_arr[ipop],0,0);
  else if(tag==BG_BZ_SRCS)
    return f_of_r_linear(par,r,par->srcs_bz_arr[ipop],par->srcs_bz_arr[ipop][0],1);
  else if(tag==BG_NORM_SRCS) {
    return f_of_r_linear(par,r,par->srcs_norm_arr[ipop],par->norm_srcs_0[ipop],
			 par->norm_srcs_f[ipop]);
  }
  else if(tag==BG_TZ_IMAP)
    return f_of_r_linear(par,r,par->imap_tz_arr[ipop],0,0);
  else if(tag==BG_BZ_IMAP)
    return f_of_r_linear(par,r,par->imap_bz_arr[ipop],par->imap_bz_arr[ipop][0],1);
  else if(tag==BG_NORM_IMAP) {
    return f_of_r_linear(par,r,par->imap_norm_arr[ipop],par->norm_imap_0[ipop],
			 par->norm_imap_f[ipop]);
  }
  else {
    report_error(1,"Unknown background quantity tag %d\n",tag);
    return -1;
  }
}

static double a_of_r_provisional(ParamCoLoRe *par,double r)
{
  if(r<=0) return 1.;
  else if(r>=par->r_arr_a2r[0]) return 1E-6;
  else {
    int ia=0;
    while(r<=par->r_arr_a2r[ia])
      ia++;
    return  par->a_arr_a2r[ia-1]+(par->a_arr_a2r[ia]-par->a_arr_a2r[ia-1])*
      (r-par->r_arr_a2r[ia-1])/(par->r_arr_a2r[ia]-par->r_arr_a2r[ia-1]);
  }
}

double r_of_z(ParamCoLoRe *par,double z)
{
  double a=1./(1+z);
  if(a>=1) return 0;
  else if(a<=0) return par->r_arr_a2r[0];
  else {
    int ia=(int)(a*(NA-1));
    double r=par->r_arr_a2r[ia]+(par->r_arr_a2r[ia+1]-par->r_arr_a2r[ia])*
      (a-par->a_arr_a2r[ia])*(NA-1.);
    return r;
  }
}

typedef struct {
  ParamCoLoRe *par;
  int l;
  double chi_0_a;
  double chi_f_a;
  double chi_0_b;
  double chi_f_b;
  double chi_s;
  int cl_type;
} ClParams;

static double window_kappa_limber(ParamCoLoRe *par,int l,double k,double chi_0,double chi_f,double chi_s)
{
  double chi_l=(l+0.5)/k;

  if((chi_l<=0) || (chi_l<chi_0) || (chi_l>chi_f))
    return 0;
  else {
    double gf=get_bg(par,chi_l,BG_D1,0);
    double aa=1/(1+get_bg(par,chi_l,BG_Z,0));
  
    return par->prefac_lensing*l*(l+1.)*gf*(chi_s-chi_l)/(k*k*chi_s*chi_l*aa);
  }
}

static double window_isw_limber(ParamCoLoRe *par,int l,double k,double chi_0,double chi_f,double chi_s)
{
  double chi_l=(l+0.5)/k;

  if((chi_l<=0) || (chi_l<chi_0) || (chi_l>chi_f))
    return 0;
  else
    return -2*par->prefac_lensing*get_bg(par,chi_l,BG_PD,0)/(k*k);
}

static double cl_integrand(double lk,void *params)
{
  ClParams *p=(ClParams *)params;
  double k=pow(10.,lk);
  double pk=pk_linear0(p->par,lk);
  double wa=0,wb=0;
  if(p->cl_type==CL_TYPE_KAPPA) {
    wa=window_kappa_limber(p->par,p->l,k,p->chi_0_a,p->chi_f_a,p->chi_s);
    wb=window_kappa_limber(p->par,p->l,k,p->chi_0_b,p->chi_f_b,p->chi_s);
  }
  else if(p->cl_type==CL_TYPE_ISW) {
    wa=window_isw_limber(p->par,p->l,k,p->chi_0_a,p->chi_f_a,p->chi_s);
    wb=window_isw_limber(p->par,p->l,k,p->chi_0_b,p->chi_f_b,p->chi_s);
  }
  if((p->par->do_smoothing) && (p->par->smooth_potential))
    pk*=exp(-k*k*p->par->r2_smooth);

  return k*wa*wb*pk;
}

static double compute_cl(ParamCoLoRe *par,int l,
			 double chi_0_a,double chi_f_a,
			 double chi_0_b,double chi_f_b,
			 double chi_s,int cl_type)
{
  ClParams p;
  gsl_function F;
  double lkmin,lkmax,eresult,result=0,chi_min=1E16,chi_max=-1;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  p.par=par;
  p.l=l;
  p.chi_0_a=chi_0_a;
  p.chi_f_a=chi_f_a;
  p.chi_0_b=chi_0_b;
  p.chi_f_b=chi_f_b;
  p.chi_s=chi_s;
  p.cl_type=cl_type;
  F.function=&cl_integrand;
  F.params=&p;

  if(chi_0_a<chi_min) chi_min=chi_0_a;
  if(chi_0_b<chi_min) chi_min=chi_0_b;
  if(chi_f_a<chi_min) chi_min=chi_f_a;
  if(chi_f_b<chi_min) chi_min=chi_f_b;
  if(chi_0_a>chi_max) chi_max=chi_0_a;
  if(chi_0_b>chi_max) chi_max=chi_0_b;
  if(chi_f_a>chi_max) chi_max=chi_f_a;
  if(chi_f_b>chi_max) chi_max=chi_f_b;
  if(chi_min<=0) lkmax=2.;
  else lkmax=log10((l+0.5)/chi_min);
  if(chi_max<=0) lkmin=-5.;
  else lkmin=log10((l+0.5)/chi_max);

  gsl_integration_qag(&F,lkmin,lkmax,0,1E-4,1000,GSL_INTEG_GAUSS41,w,&result,&eresult);
  gsl_integration_workspace_free(w);

  return M_LN10*result/(l+0.5);
}

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

static void compute_hp_shell_distances_imap(HealpixShells *shell,flouble nu_rest,
					    char *fname_nutable,Csm_params *pars)
{
  int ii;
  FILE *fi;

  //Figure out radial shells
  fi=fopen(fname_nutable,"r");
  if(fi==NULL) error_open_file(fname_nutable);
  for(ii=0;ii<shell->nr;ii++) {
    double nu0,nuf;
    int stat=fscanf(fi,"%lf %lf",&nu0,&nuf);
    if(stat!=2) error_read_line(fname_nutable,ii+1);
    shell->r0[ii]=csm_radial_comoving_distance(pars,nuf/nu_rest);
    shell->rf[ii]=csm_radial_comoving_distance(pars,nu0/nu_rest);
  }
  fclose(fi);
}

void cosmo_set(ParamCoLoRe *par)
{
  //////
  // This initializes the cosmological model
  // at redshift z_s

  int ii,ipop,nz;
  double *zarr,*fzarr;
  FILE *fi;
  gsl_spline *spline_srcs_bz[NPOP_MAX];
  gsl_spline *spline_srcs_nz[NPOP_MAX];
  gsl_interp_accel *intacc_srcs=gsl_interp_accel_alloc();
  gsl_spline *spline_imap_bz[NPOP_MAX];
  gsl_spline *spline_imap_tz[NPOP_MAX];
  gsl_interp_accel *intacc_imap=gsl_interp_accel_alloc();

  Csm_params *pars=csm_params_new();
  csm_unset_gsl_eh();
  csm_background_set(pars,par->OmegaM,par->OmegaL,par->OmegaB,par->weos,0,par->hhub,2.275);

  par->fgrowth_0=csm_f_growth(pars,1);
  par->hubble_0=csm_hubble(pars,1);

  par->r_min=csm_radial_comoving_distance(pars,1/(1+par->z_min));
  par->r_max=csm_radial_comoving_distance(pars,1/(1+par->z_max));

  par->prefac_lensing=1.5*par->hubble_0*par->hubble_0*par->OmegaM;

  par->l_box=2*par->r_max*(1+2./par->n_grid);
  par->pos_obs[0]=0.5*par->l_box;
  par->pos_obs[1]=0.5*par->l_box;
  par->pos_obs[2]=0.5*par->l_box;

  for(ipop=0;ipop<par->n_srcs;ipop++) {
    fi=fopen(par->fnameBzSrcs[ipop],"r");
    if(fi==NULL) error_open_file(par->fnameBzSrcs[ipop]);
    nz=linecount(fi); rewind(fi);
    zarr=my_malloc(nz*sizeof(double));
    fzarr=my_malloc(nz*sizeof(double));
    for(ii=0;ii<nz;ii++) {
      int stat=fscanf(fi,"%lf %lf",&(zarr[ii]),&(fzarr[ii]));
      if(stat!=2) error_read_line(par->fnameBzSrcs[ipop],ii+1);
    }
    if((zarr[0]>par->z_min) || (zarr[nz-1]<par->z_max))
      report_error(1,"Bias z-range is too small\n");
    spline_srcs_bz[ipop]=gsl_spline_alloc(gsl_interp_cspline,nz);
    gsl_spline_init(spline_srcs_bz[ipop],zarr,fzarr,nz);
    free(zarr); free(fzarr);
    fclose(fi);

    fi=fopen(par->fnameNzSrcs[ipop],"r");
    if(fi==NULL) error_open_file(par->fnameNzSrcs[ipop]);
    nz=linecount(fi); rewind(fi);
    zarr=my_malloc(nz*sizeof(double));
    fzarr=my_malloc(nz*sizeof(double));
    for(ii=0;ii<nz;ii++) {
      double rz,hz,a;
      int stat=fscanf(fi,"%lf %lf",&(zarr[ii]),&(fzarr[ii]));
      if(stat!=2) error_read_line(par->fnameNzSrcs[ipop],ii+1);
      a=1./(1+zarr[ii]);
      hz=csm_hubble(pars,a);
      rz=csm_radial_comoving_distance(pars,a);
      fzarr[ii]*=RTOD*RTOD*hz/(rz*rz);
    }
    //Correct for z[0]=0
    if(zarr[0]==0)
      fzarr[0]=fzarr[1];
    if((zarr[0]>par->z_min) || (zarr[nz-1]<par->z_max))
      report_error(1,"N(z) z-range is too small\n");
    spline_srcs_nz[ipop]=gsl_spline_alloc(gsl_interp_cspline,nz);
    gsl_spline_init(spline_srcs_nz[ipop],zarr,fzarr,nz);
    free(zarr); free(fzarr);
    fclose(fi);

    par->srcs_nz_arr[ipop]=my_malloc(NA*sizeof(double));
    par->srcs_bz_arr[ipop]=my_malloc(NA*sizeof(double));
  }

  for(ipop=0;ipop<par->n_imap;ipop++) {
    fi=fopen(par->fnameBzImap[ipop],"r");
    if(fi==NULL) error_open_file(par->fnameBzImap[ipop]);
    nz=linecount(fi); rewind(fi);
    zarr=my_malloc(nz*sizeof(double));
    fzarr=my_malloc(nz*sizeof(double));
    for(ii=0;ii<nz;ii++) {
      int stat=fscanf(fi,"%lf %lf",&(zarr[ii]),&(fzarr[ii]));
      if(stat!=2) error_read_line(par->fnameBzImap[ipop],ii+1);
    }
    if((zarr[0]>par->z_min) || (zarr[nz-1]<par->z_max))
      report_error(1,"Bias z-range is too small\n");
    spline_imap_bz[ipop]=gsl_spline_alloc(gsl_interp_cspline,nz);
    gsl_spline_init(spline_imap_bz[ipop],zarr,fzarr,nz);
    free(zarr); free(fzarr);
    fclose(fi);

    fi=fopen(par->fnameTzImap[ipop],"r");
    if(fi==NULL) error_open_file(par->fnameTzImap[ipop]);
    nz=linecount(fi); rewind(fi);
    zarr=my_malloc(nz*sizeof(double));
    fzarr=my_malloc(nz*sizeof(double));
    for(ii=0;ii<nz;ii++) {
      int stat=fscanf(fi,"%lf %lf",&(zarr[ii]),&(fzarr[ii]));
      if(stat!=2) error_read_line(par->fnameTzImap[ipop],ii+1);
    }
    if((zarr[0]>par->z_min) || (zarr[nz-1]<par->z_max))
      report_error(1,"T(z) z-range is too small\n");
    spline_imap_tz[ipop]=gsl_spline_alloc(gsl_interp_cspline,nz);
    gsl_spline_init(spline_imap_tz[ipop],zarr,fzarr,nz);
    free(zarr); free(fzarr);
    fclose(fi);

    par->imap_tz_arr[ipop]=my_malloc(NA*sizeof(double));
    par->imap_bz_arr[ipop]=my_malloc(NA*sizeof(double));
  }

  //Set z-dependent functions
  for(ii=0;ii<NA;ii++) {
    double a=((double)ii)/(NA-1);
    double ra=csm_radial_comoving_distance(pars,a);
    par->a_arr_a2r[ii]=a;
    par->r_arr_a2r[ii]=ra;
  }

  double growth0=csm_growth_factor(pars,1);
  par->glob_idr=(NA-1)/par->r_arr_a2r[0];
  for(ii=0;ii<NA;ii++) {
    double r=ii/par->glob_idr;
    double a=a_of_r_provisional(par,r);
    double z=1./a-1;
    double gz=csm_growth_factor(pars,a)/growth0;
    double om=csm_omega_m(pars,a);
    double fz=csm_f_growth(pars,a);
    double hhz=csm_hubble(pars,a);
    par->z_arr_r2z[ii]=z;
    par->r_arr_r2z[ii]=r;
    par->growth_d_arr[ii]=gz;
    par->growth_d2_arr[ii]=-0.42857142857*gz*gz*pow(om,-0.00699300699);
    par->growth_v_arr[ii]=(gz*hhz*fz)/(par->fgrowth_0*par->hubble_0); //This is for the comoving velocity
    par->growth_pd_arr[ii]=gz*hhz*(fz-1);
    par->ihub_arr[ii]=1./hhz;
    for(ipop=0;ipop<par->n_srcs;ipop++) {
      if((z<par->z_min) || (z>par->z_max)) {
	par->srcs_nz_arr[ipop][ii]=0;
	par->srcs_bz_arr[ipop][ii]=1.;
      }
      par->srcs_nz_arr[ipop][ii]=gsl_spline_eval(spline_srcs_nz[ipop],z,intacc_srcs);
      par->srcs_bz_arr[ipop][ii]=gsl_spline_eval(spline_srcs_bz[ipop],z,intacc_srcs);
    }
    for(ipop=0;ipop<par->n_imap;ipop++) {
      if((z<par->z_min) || (z>par->z_max)) {
	par->imap_tz_arr[ipop][ii]=0;
	par->imap_bz_arr[ipop][ii]=1.;
      }
      par->imap_tz_arr[ipop][ii]=gsl_spline_eval(spline_imap_tz[ipop],z,intacc_imap);
      par->imap_bz_arr[ipop][ii]=gsl_spline_eval(spline_imap_bz[ipop],z,intacc_imap);
    }
  }

  gsl_interp_accel_free(intacc_srcs);
  for(ipop=0;ipop<par->n_imap;ipop++) {
    gsl_spline_free(spline_imap_tz[ipop]);
    gsl_spline_free(spline_imap_bz[ipop]);
  }
  gsl_interp_accel_free(intacc_imap);
  for(ipop=0;ipop<par->n_srcs;ipop++) {
    gsl_spline_free(spline_srcs_nz[ipop]);
    gsl_spline_free(spline_srcs_bz[ipop]);
  }

  pk_linear_set(par);

  if(par->do_kappa) {
    for(ii=0;ii<par->n_kappa;ii++) {
      double z=par->z_kappa_out[ii];
      if(z>par->z_max) {
#ifdef _ADD_EXTRA_KAPPA
	int l,nl=3*par->nside_kappa;
	double chi_here=r_of_z(par,z);
#ifdef _DEBUG
	if(NodeThis==0)
	  fprintf(par->f_dbg,"Power spectra for extra kappa, shell %d (z=%.3lf)\n",ii+1,z);
#endif //_DEBUG
	par->need_extra_kappa[ii]=1;
	par->fl_mean_extra_kappa[ii]=my_malloc(nl*sizeof(flouble));
	par->cl_extra_kappa[ii]=my_malloc(nl*sizeof(flouble));
	for(l=0;l<nl;l++) {
	  flouble cl_k_11=compute_cl(par,l,0.        ,par->r_max,0.        ,par->r_max,chi_here,CL_TYPE_KAPPA);
	  flouble cl_k_12=compute_cl(par,l,0.        ,par->r_max,par->r_max,chi_here  ,chi_here,CL_TYPE_KAPPA);
	  flouble cl_k_22=compute_cl(par,l,par->r_max,chi_here  ,par->r_max,chi_here  ,chi_here,CL_TYPE_KAPPA);
	  par->fl_mean_extra_kappa[ii][l]=cl_k_12/cl_k_11;
	  par->cl_extra_kappa[ii][l]=cl_k_22-cl_k_12*cl_k_12/cl_k_11;
	  if(l==0) {
	    par->fl_mean_extra_kappa[ii][l]=0;
	    par->cl_extra_kappa[ii][l]=0;
	  }
#ifdef _DEBUG
	  if(NodeThis==0)
	    fprintf(par->f_dbg,"%d %lE %lE \n",l,par->fl_mean_extra_kappa[ii][l],par->cl_extra_kappa[ii][l]);
#endif //_DEBUG
	}
#ifdef _DEBUG
	if(NodeThis==0)
	  fprintf(par->f_dbg,"\n");
#endif //_DEBUG
#else //_ADD_EXTRA_KAPPA
	report_error(1,"Source plane %d outside redshift range\n",ii+1);
#endif //_ADD_EXTRA_KAPPA
      }
    }
  }

  if(par->do_isw) {
    for(ii=0;ii<par->n_isw;ii++) {
      double z=par->z_isw_out[ii];
      if(z>par->z_max) {
#ifdef _ADD_EXTRA_KAPPA
	int l,nl=3*par->nside_isw;
	double chi_here=r_of_z(par,z);
#ifdef _DEBUG
	if(NodeThis==0)
	  fprintf(par->f_dbg,"Power spectra for extra isw, shell %d (z=%.3lf)\n",ii+1,z);
#endif //_DEBUG
	par->need_extra_isw[ii]=1;
	par->fl_mean_extra_isw[ii]=my_malloc(nl*sizeof(flouble));
	par->cl_extra_isw[ii]=my_malloc(nl*sizeof(flouble));
	for(l=0;l<nl;l++) {
	  flouble cl_i_11=compute_cl(par,l,0.        ,par->r_max,0.        ,par->r_max,chi_here,CL_TYPE_ISW);
	  flouble cl_i_12=compute_cl(par,l,0.        ,par->r_max,par->r_max,chi_here  ,chi_here,CL_TYPE_ISW);
	  flouble cl_i_22=compute_cl(par,l,par->r_max,chi_here  ,par->r_max,chi_here  ,chi_here,CL_TYPE_ISW);
	  par->fl_mean_extra_isw[ii][l]=cl_i_12/cl_i_11;
	  par->cl_extra_isw[ii][l]=cl_i_22-cl_i_12*cl_i_12/cl_i_11;
	  if(l==0) {
	    par->fl_mean_extra_isw[ii][l]=0;
	    par->cl_extra_isw[ii][l]=0;
	  }
#ifdef _DEBUG
	  if(NodeThis==0)
	    fprintf(par->f_dbg,"%d %lE %lE \n",l,par->fl_mean_extra_isw[ii][l],par->cl_extra_isw[ii][l]);
#endif //_DEBUG
	}
#ifdef _DEBUG
	if(NodeThis==0)
	  fprintf(par->f_dbg,"\n");
#endif //_DEBUG
#else //_ADD_EXTRA_KAPPA
	report_error(1,"Source plane %d outside redshift range\n",ii+1);
#endif //_ADD_EXTRA_KAPPA
      }
    }
  }

  csm_params_free(pars);
  //FREEE SLPINES!
}

void compute_tracer_cosmo(ParamCoLoRe *par)
{
  Csm_params *pars=csm_params_new();
  csm_unset_gsl_eh();
  csm_set_verbosity(0);
  csm_background_set(pars,par->OmegaM,par->OmegaL,par->OmegaB,par->weos,0,par->hhub,2.275);

  if(par->do_imap) {
    int ipop;
    for(ipop=0;ipop<par->n_imap;ipop++) {
      compute_hp_shell_distances_imap(par->imap[ipop],par->nu0_imap[ipop],par->fnameNuImap[ipop],pars);
    }
  }
  if(par->do_kappa) {
    int ii;
    for(ii=0;ii<par->n_kappa;ii++) {
      double z=par->z_kappa_out[ii];
      par->kmap->r0[ii]=csm_radial_comoving_distance(pars,1./(1+z));
      par->kmap->rf[ii]=csm_radial_comoving_distance(pars,1./(1+z));
    }
  }
  if(par->do_isw) {
    int ii;
    for(ii=0;ii<par->n_isw;ii++) {
      double z=par->z_isw_out[ii];
      par->pd_map->r0[ii]=csm_radial_comoving_distance(pars,1./(1+z));
      par->pd_map->rf[ii]=csm_radial_comoving_distance(pars,1./(1+z));
    }
  }

  csm_params_free(pars);
}
