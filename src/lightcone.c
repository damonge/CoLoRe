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

static inline double get_cosine(double index,double dx)
{
#if PIXTYPE==PT_CEA
  return index*dx-1;
#elif PIXTYPE==PT_CAR
  return cos(M_PI-index*dx);
#endif //PIXTYPE
}

static void get_sources_single(ParamCoLoRe *par,int ipop)
{
  //////
  // Uses the gaussian matter density field to obtain a
  // poisson sampling of point sources (returns an integer array
  // with the number of sources in each cell).

  int ii,nthr;
  lint *np_tot_thr;
#ifdef _HAVE_OMP
  nthr=omp_get_max_threads();
#else //_HAVE_OMP
  nthr=1;
#endif //_HAVE_OMP
  np_tot_thr=my_calloc(nthr,sizeof(lint));

  print_info(" %d-th galaxy population\n",ipop);
  if(NodeThis==0) timer(0);
  print_info("   Poisson-sampling\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,np_tot_thr,IThread0,ipop,nthr)
#endif //_HAVE_OMP
  {
    int ir;
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_srcs*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);

#ifdef _HAVE_OMP
#pragma omp for schedule(static) //TODO: this will give very bad load balance
#endif //_HAVE_OMP 
    for(ir=0;ir<par->oi_beams[0]->nr;ir++) {
      double r0=par->oi_beams[0]->r0_arr[ir];
      double rf=par->oi_beams[0]->rf_arr[ir];
      double rm=(rf+r0)*0.5;
      double redshift=z_of_r(par,rm);
      double ndens=ndens_of_z_srcs(par,redshift,ipop);
      if(ndens>0) {
	int ib;
	int nside=par->oi_beams[0]->nside_arr[ir];
	double dphi=M_PI/nside;
#if PIXTYPE==PT_CEA
	double dcth=2./nside;
	double cell_vol=(rf*rf*rf-r0*r0*r0)*dcth*dphi/3;
#elif PIXTYPE==PT_CAR
	double dth=M_PI/nside;
#endif //PIXTYPE
	double bias=bias_of_z_srcs(par,redshift,ipop);
	double gfb=dgrowth_of_r(par,rm)*bias;
	for(ib=0;ib<par->n_beams_here;ib++) {
	  int ind_cth;
	  int ncth=par->oi_beams[ib]->icthf_arr[ir]-par->oi_beams[ib]->icth0_arr[ir]+1;
	  int nphi=par->oi_beams[ib]->iphif_arr[ir]-par->oi_beams[ib]->iphi0_arr[ir]+1;
	  flouble *dens_slice=par->dens_beams[ib][ir];
	  int *nsrc_slice=par->nsrc_beams[ib][ir];
	  for(ind_cth=0;ind_cth<ncth;ind_cth++) {
	    int ind_phi;
	    int ind_cth_t=ind_cth*nphi;
#if PIXTYPE==PT_CAR
	    double dcth=get_cosine(par->oi_beams[ib]->icth0_arr[ir]+ind_cth+1.0,dth)-
	      get_cosine(par->oi_beams[ib]->icth0_arr[ir]+ind_cth+0.0,dth);
	    double cell_vol=(rf*rf*rf-r0*r0*r0)*dcth*dphi/3;
#endif //PIXTYPE
	    for(ind_phi=0;ind_phi<nphi;ind_phi++) {
	      double lambda=ndens*cell_vol*exp(gfb*(dens_slice[ind_cth_t+ind_phi]-0.5*gfb*par->sigma2_gauss));
	      int npp=rng_poisson(lambda,rng_thr);
	      nsrc_slice[ind_cth_t+ind_phi]=npp;
	      np_tot_thr[ithr]+=npp;
	    }
	  }
	}
      }
    } //end omp for
    end_rng(rng_thr);
  } //end omp parallel

  par->nsources_this[ipop]=0;
  for(ii=0;ii<nthr;ii++)
    par->nsources_this[ipop]+=np_tot_thr[ii];

  lint nsources_total=0;
#ifdef _HAVE_MPI
  MPI_Allreduce(&(par->nsources_this[ipop]),&nsources_total,1,LINT_MPI,MPI_SUM,MPI_COMM_WORLD);
#else //_HAVE_MPI
  nsources_total=par->nsources_this[ipop];
#endif //_HAVE_MPI

  print_info("   There will be %ld objects in total \n",(long)nsources_total);
#ifdef _DEBUG
  fprintf(par->f_dbg,"Node %d has %ld particles\n",NodeThis,(long)(par->nsources_this));
#endif //_DEBUG

  for(ii=nthr-1;ii>0;ii--) {
    int jj;
    lint nh=0;
    for(jj=0;jj<ii;jj++)
      nh+=np_tot_thr[jj];
    np_tot_thr[ii]=nh;
  }
  np_tot_thr[0]=0;
  //np_tot_thr now contains the id of the first particle in the thread
  
  par->srcs[ipop]=my_malloc(par->nsources_this[ipop]*sizeof(Src));

  if(NodeThis==0) timer(0);
  print_info("   Assigning coordinates\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none)			\
  shared(par,IThread0,np_tot_thr,ipop,nthr)
#endif //_HAVE_OMP
  {
    int ir;
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_srcs*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);

#ifdef _HAVE_OMP
#pragma omp for schedule(static) //TODO: this will give very bad load balance
#endif //_HAVE_OMP
    for(ir=0;ir<par->oi_beams[0]->nr;ir++) {
      double r0=par->oi_beams[0]->r0_arr[ir];
      double rf=par->oi_beams[0]->rf_arr[ir];
      double dr=rf-r0;
      double rm=r0+0.5*dr;
      double redshift=z_of_r(par,rm);
      double ndens=ndens_of_z_srcs(par,redshift,ipop);
      if(ndens>0) {
	int ib;
	int nside=par->oi_beams[0]->nside_arr[ir];
#if PIXTYPE==PT_CEA
	double dcth=2./nside;
#elif PIXTYPE==PT_CAR
	double dth=M_PI/nside;
#endif //PIXTYPE
	double dphi=M_PI/nside;
	double vg=vgrowth_of_r(par,rm);
	for(ib=0;ib<par->n_beams_here;ib++) {
	  int ind_cth;
	  int icth0=par->oi_beams[ib]->icth0_arr[ir];
	  int iphi0=par->oi_beams[ib]->iphi0_arr[ir];
	  int ncth=par->oi_beams[ib]->icthf_arr[ir]-par->oi_beams[ib]->icth0_arr[ir]+1;
	  int nphi=par->oi_beams[ib]->iphif_arr[ir]-par->oi_beams[ib]->iphi0_arr[ir]+1;
	  flouble *vrad_slice=par->vrad_beams[ib][ir];
	  int *nsrc_slice=par->nsrc_beams[ib][ir];
	  for(ind_cth=0;ind_cth<ncth;ind_cth++) {
	    int ind_phi;
	    int ind_cth_t=ind_cth*nphi;
#if PIXTYPE==PT_CEA
	    double cth0=get_cosine(ind_cth+icth0+0.0,dcth);
#elif PIXTYPE==PT_CAR
	    double cth0=get_cosine(ind_cth+icth0+0.0,dth);
	    double dcth=get_cosine(ind_cth+icth0+1.0,dth)-cth0;
#endif //PIXTYPE
	    for(ind_phi=0;ind_phi<nphi;ind_phi++) {
	      int ip;
	      double e1=0,e2=0;
	      double phi0=(ind_phi+iphi0)*dphi;
	      int npp=nsrc_slice[ind_cth_t+ind_phi];
	      double dz_rsd=vg*vrad_slice[ind_cth_t+ind_phi];
	      if(par->shear_srcs[ipop]) {
		double pxx=par->p_xx_beams[ib][ir][ind_cth_t+ind_phi];
		double pxy=par->p_xy_beams[ib][ir][ind_cth_t+ind_phi];
		double pyy=par->p_yy_beams[ib][ir][ind_cth_t+ind_phi];
		double g1=pxx-pyy;
		double g2=2*pxy;
#ifdef _NONLINEAR_ELLIPTICITIES
		double kappa=pxx+pyy;
		double fac=2*(1-kappa)/((1-kappa)*(1-kappa)+g1*g1+g2*g2);
#else //_NONLINEAR_ELLIPTICITIES
		double fac=2;
#endif //_NONLINEAR_ELLIPTICITIES
		e1=fac*g1;
		e2=fac*g2;
	      }
	      for(ip=0;ip<npp;ip++) {
		double cth=cth0+dcth*rng_01(rng_thr);
		double phi=phi0+dphi*rng_01(rng_thr);
		double r=r0+dr*rng_01(rng_thr); //TODO: this is not completely correct
		lint pid=np_tot_thr[ithr];
		par->srcs[ipop][pid].ra=RTOD*phi;
		par->srcs[ipop][pid].dec=90-RTOD*acos(cth);
		par->srcs[ipop][pid].z0=z_of_r(par,r);
		par->srcs[ipop][pid].dz_rsd=dz_rsd;
		par->srcs[ipop][pid].e1=e1;
		par->srcs[ipop][pid].e2=e2;
		np_tot_thr[ithr]++;
	      }
	    }
	  }
	}
      }
    } //end omp for
    end_rng(rng_thr);
  } //end omp parallel
  if(NodeThis==0) timer(2);

  free(np_tot_thr);
}

//////
// Integrates the Newtonian potential along the line of sight
// with the lensing kernel.
void integrate_lensing(ParamCoLoRe *par)
{
  int ib;

#ifdef _DEBUG
  print_info("*** Integrating lensing\n");
  if(NodeThis==0) timer(0);
#endif //_DEBUG
  for(ib=0;ib<par->n_beams_here;ib++) {
    int ir;
    int nside_old=par->oi_beams[ib]->nside_arr[0];
    double *pxx1_old=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pxx1_new=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pxx2_old=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pxx2_new=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pxy1_old=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pxy1_new=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pxy2_old=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pxy2_new=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pyy1_old=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pyy1_new=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pyy2_old=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    double *pyy2_new=my_calloc(par->oi_beams[ib]->num_pix[par->oi_beams[ib]->nr-1],sizeof(double));
    for(ir=0;ir<par->oi_beams[ib]->nr;ir++) {
      int icth;
      int nside_new=par->oi_beams[ib]->nside_arr[ir];
      int nside_ratio=nside_new/nside_old;
      int ncth=par->oi_beams[ib]->icthf_arr[ir]-par->oi_beams[ib]->icth0_arr[ir]+1;
      int nphi=par->oi_beams[ib]->iphif_arr[ir]-par->oi_beams[ib]->iphi0_arr[ir]+1;
      int nphi_old=nphi/nside_ratio;
      double r0=par->oi_beams[ib]->r0_arr[ir];
      double rf=par->oi_beams[ib]->rf_arr[ir];
      double rm=0.5*(r0+rf),dr=rf-r0;
      double redshift=z_of_r(par,rm);
      double g_phi=dgrowth_of_r(par,rm)*(1+redshift);
      //      double integ1=g_phi*0.5*(rf*rf-r0*r0);
      //      double integ2=g_phi*(rf*rf*rf-r0*r0*r0)/3;
      double integ1=g_phi*rm*dr;
      double integ2=g_phi*rm*rm*dr;
      if(ncth*nphi!=par->oi_beams[ib]->num_pix[ir])
	report_error(1,"WTF\n");
      for(icth=0;icth<ncth;icth++) {
	int iphi;
	int icth_old=icth/nside_ratio;
	for(iphi=0;iphi<nphi;iphi++) {
	  int iphi_old=iphi/nside_ratio;
	  int index_new=iphi+nphi*icth;
	  int index_old=iphi_old+nphi_old*icth_old;
	  pxx1_new[index_new]=pxx1_old[index_old]+integ1*par->p_xx_beams[ib][ir][index_new];
	  pxx2_new[index_new]=pxx2_old[index_old]+integ2*par->p_xx_beams[ib][ir][index_new];
	  par->p_xx_beams[ib][ir][index_new]=pxx1_new[index_new]-pxx2_new[index_new]/rf;
	  pxy1_new[index_new]=pxy1_old[index_old]+integ1*par->p_xy_beams[ib][ir][index_new];
	  pxy2_new[index_new]=pxy2_old[index_old]+integ2*par->p_xy_beams[ib][ir][index_new];
	  par->p_xy_beams[ib][ir][index_new]=pxy1_new[index_new]-pxy2_new[index_new]/rf;
	  pyy1_new[index_new]=pyy1_old[index_old]+integ1*par->p_yy_beams[ib][ir][index_new];
	  pyy2_new[index_new]=pyy2_old[index_old]+integ2*par->p_yy_beams[ib][ir][index_new];
	  par->p_yy_beams[ib][ir][index_new]=pyy1_new[index_new]-pyy2_new[index_new]/rf;
	}
      }
      nside_old=nside_new;
      memcpy(pxx1_old,pxx1_new,par->oi_beams[ib]->num_pix[ir]*sizeof(double));
      memcpy(pxx2_old,pxx2_new,par->oi_beams[ib]->num_pix[ir]*sizeof(double));
      memcpy(pxy1_old,pxy1_new,par->oi_beams[ib]->num_pix[ir]*sizeof(double));
      memcpy(pxy2_old,pxy2_new,par->oi_beams[ib]->num_pix[ir]*sizeof(double));
      memcpy(pyy1_old,pyy1_new,par->oi_beams[ib]->num_pix[ir]*sizeof(double));
      memcpy(pyy2_old,pyy2_new,par->oi_beams[ib]->num_pix[ir]*sizeof(double));
    }
    free(pxx1_new); free(pxx1_old);
    free(pxx2_new); free(pxx2_old);
    free(pxy1_new); free(pxy1_old);
    free(pxy2_new); free(pxy2_old);
    free(pyy1_new); free(pyy1_old);
    free(pyy2_new); free(pyy2_old);
  }
#ifdef _DEBUG
  if(NodeThis==0) timer(2);
#endif //_DEBUG
}

void get_sources(ParamCoLoRe *par)
{
  int ipop;

  //First, compute lensing Hessian
  print_info("*** Getting point sources\n");
  for(ipop=0;ipop<par->n_srcs;ipop++)
    get_sources_single(par,ipop);
  print_info("\n");
}

static int get_r_index_imap(HealpixShells *sh,double r,int ir_start)
{
  int gotit=0;
  int ir0;
  if(ir_start<0)
    ir0=0;
  else if(ir_start>=sh->nr)
    ir0=sh->nr-1;
  else
    ir0=ir_start;

  while(!gotit) {
    if((ir0==-1) || (ir0==sh->nr))
      gotit=1;
    else {
      if(r<sh->r0[ir0])
	ir0++;
      else {
        if(r>=sh->rf[ir0])
          ir0--;
	else
          gotit=1;
      }
    }
  }

  return ir0;
}

static void find_shell_pixels(ParamCoLoRe *par,HealpixShells *shell)
{
  long ip,npx=he_nside2npix(shell->nside);
  double pixsize=2*sqrt(4*M_PI/npx);
  shell->num_pix=0;
  shell->listpix=my_malloc(npx*sizeof(long));

  for(ip=0;ip<npx;ip++) {
    int ib;
    int goodpix=0;
    double phi0,phif,phim;
    double th0,thf,thm,cth0,cthf;
    pix2ang_ring(shell->nside,ip,&thm,&phim);
    phi0=phim-pixsize;
    phif=phim+pixsize;
    th0=CLAMP(thm-pixsize,0,M_PI);
    thf=CLAMP(thm+pixsize,0,M_PI);
    cthf=cos(th0);
    cth0=cos(thf);
    for(ib=0;ib<par->n_beams_here;ib++) {
      OnionInfo *beam=par->oi_beams[ib];
      double cth0_b,cthf_b,phi0_b,phif_b;
#if PIXTYPE==PT_CEA
      cth0_b=get_cosine(beam->icth0_arr[0]+0.0,2./beam->nside_arr[0]);
      cthf_b=get_cosine(beam->icth0_arr[0]+1.0,2./beam->nside_arr[0]);
#elif PIXTYPE==PT_CAR
      cth0_b=get_cosine(beam->icth0_arr[0]+0.0,M_PI/beam->nside_arr[0]);
      cthf_b=get_cosine(beam->icth0_arr[0]+1.0,M_PI/beam->nside_arr[0]);
#endif //PIXTYPE
      phi0_b=beam->iphi0_arr[0]*M_PI/beam->nside_arr[0];
      phif_b=(beam->iphif_arr[0]+1)*M_PI/beam->nside_arr[0];
      if(((cth0<=cthf_b)&&(cth0>=cth0_b)) || ((cthf<=cthf_b)&&(cthf>=cth0_b))) { //cth in range
	if(((phi0<=phif_b)&&(phi0>=phi0_b)) || ((phif<=phif_b)&&(phif>=phi0_b))) //phi in range
	  goodpix=1;
      }
    }
    if(goodpix) {
      shell->listpix[ip]=shell->num_pix;
      shell->num_pix++;
    }
    else
      shell->listpix[ip]=-1;
  }

  shell->data=my_calloc(shell->nr*shell->num_pix,sizeof(flouble));
}

#define N_SUBVOL 10
static void get_imap_single(ParamCoLoRe *par,int ipop)
{
  int nthr;
#ifdef _HAVE_OMP
  nthr=omp_get_max_threads();
#else //_HAVE_OMP
  nthr=1;
#endif //_HAVE_OMP

  print_info(" %d-th IM species\n",ipop);
  if(NodeThis==0) timer(0);
  find_shell_pixels(par,par->imap[ipop]);

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,IThread0,ipop,nthr)
#endif //_HAVE_OMP
  {
    int ir;
    double hpix_area=4*M_PI/he_nside2npix(par->imap[ipop]->nside);
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_imap*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);
    double *disp=my_calloc(3*N_SUBVOL,sizeof(double));
    for(ir=0;ir<3*N_SUBVOL;ir++)
      disp[ir]=rng_01(rng_thr);
    end_rng(rng_thr);

#ifdef _HAVE_OMP
#pragma omp for schedule(dynamic)
#endif //_HAVE_OMP 
    for(ir=0;ir<par->oi_beams[0]->nr;ir++) {
      double r0=par->oi_beams[0]->r0_arr[ir];
      double rf=par->oi_beams[0]->rf_arr[ir];
      double rm=(rf+r0)*0.5;
      double dr=rf-r0;
      double redshift=z_of_r(par,rm);
      double tmean=temp_of_z_imap(par,redshift,ipop);
      if(tmean>0) {
	int ib;
	int irad=0;
	int nside=par->oi_beams[0]->nside_arr[ir];
	double dphi=M_PI/nside;
#if PIXTYPE==PT_CEA
	double dcth=2./nside;
	double cell_vol=(rf*rf*rf-r0*r0*r0)*dcth*dphi/3/N_SUBVOL;
#elif PIXTYPE==PT_CAR
	double dth=M_PI/nside;
#endif //PIXTYPE
	double bias=bias_of_z_imap(par,redshift,ipop);
	double gfb=dgrowth_of_r(par,rm)*bias;
	double prefac_rsd=ihub_of_r(par,rm)*vgrowth_of_r(par,rm);
	for(ib=0;ib<par->n_beams_here;ib++) {
	  int ind_cth;
	  int icth0=par->oi_beams[ib]->icth0_arr[ir];
	  int iphi0=par->oi_beams[ib]->iphi0_arr[ir];
	  int ncth=par->oi_beams[ib]->icthf_arr[ir]-par->oi_beams[ib]->icth0_arr[ir]+1;
	  int nphi=par->oi_beams[ib]->iphif_arr[ir]-par->oi_beams[ib]->iphi0_arr[ir]+1;
	  flouble *dens_slice=par->dens_beams[ib][ir];
	  flouble *vrad_slice=par->vrad_beams[ib][ir];
	  for(ind_cth=0;ind_cth<ncth;ind_cth++) {
	    int ind_phi;
	    int ind_cth_t=ind_cth*nphi;
#if PIXTYPE==PT_CEA
	    double cth0=get_cosine(ind_cth+icth0+0.0,dcth);
#elif PIXTYPE==PT_CAR
	    double cth0=get_cosine(ind_cth+icth0+0.0,dth);
	    double dcth=get_cosine(ind_cth+icth0+1.0,dth)-cth0;
	    double cell_vol=(rf*rf*rf-r0*r0*r0)*dcth*dphi/3/N_SUBVOL;
#endif //PIXTYPE
	    for(ind_phi=0;ind_phi<nphi;ind_phi++) {
	      int ip;
	      double phi0=(ind_phi+iphi0)*dphi;
	      double temp=tmean*cell_vol*exp(gfb*(dens_slice[ind_cth_t+ind_phi]-0.5*gfb*par->sigma2_gauss));
	      double dr_rsd=prefac_rsd*vrad_slice[ind_cth_t+ind_phi];
	      for(ip=0;ip<N_SUBVOL;ip++) {
		double r=r0+dr_rsd+dr*disp[3*ip+0];
		irad=get_r_index_imap(par->imap[ipop],r,irad);
		if((irad>=0) && (irad<par->imap[ipop]->nr)) {
		  long pix_id;
		  long irad_t=irad*par->imap[ipop]->num_pix;
		  double cth=cth0+dcth*disp[3*ip+1];
		  double phi=phi0+dphi*disp[3*ip+2];
		  long ipix=he_ang2pix(par->imap[ipop]->nside,cth,phi);
		  pix_id=par->imap[ipop]->listpix[ipix];
		  if(pix_id<0)
		    report_error(1,"NOOOO\n");
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
		  par->imap[ipop]->data[irad_t+pix_id]+=temp;
		}
	      }
	    }
	  }
	}
      }
    } //end omp for
    free(disp);

#ifdef _HAVE_OMP
#pragma omp for schedule(dynamic)
#endif //_HAVE_OMP 
    for(ir=0;ir<par->imap[ipop]->nr;ir++) {
      long ipix;
      double r0=par->imap[ipop]->r0[ir];
      double rf=par->imap[ipop]->rf[ir];
      double i_pixel_vol=1./((rf*rf*rf-r0*r0*r0)*hpix_area/3);
      long ir_t=ir*par->imap[ipop]->num_pix;
      for(ipix=0;ipix<par->imap[ipop]->num_pix;ipix++) {
	long index=ir_t+ipix;
	par->imap[ipop]->data[index]*=i_pixel_vol;
      }
    }//end omp for
  } //end omp parallel

  if(NodeThis==0) timer(2);
}

void get_imap(ParamCoLoRe *par)
{
  int ipop;

  //First, compute lensing Hessian
  print_info("*** Getting intensity maps\n");
  for(ipop=0;ipop<par->n_imap;ipop++)
    get_imap_single(par,ipop);
  print_info("\n");
}

void get_kappa(ParamCoLoRe *par)
{

  print_info("*** Getting kappa maps\n");
  if(NodeThis==0) timer(0);
  find_shell_pixels(par,par->kmap);

#pragma omp parallel default(none) \
  shared(par)
  {
    int ir,ipx;
    double inv_hpix_area=he_nside2npix(par->kmap->nside)/(4*M_PI);
    int *nparr=malloc(par->kmap->num_pix*sizeof(int));

    //Maybe OMP this
#pragma omp for
    for(ir=0;ir<par->n_kappa;ir++) {
      int ib,irb=0;
      double r=par->kmap->r0[ir];
      long irad_t=ir*par->kmap->num_pix;
      while(irb<par->oi_beams[0]->nr) {
	if((r>=par->oi_beams[0]->r0_arr[irb]) &&
	   (r<=par->oi_beams[0]->rf_arr[irb]))
	  break;
	else
	  irb++;
      }
      if(irb>=par->oi_beams[0]->nr) {
	irb=par->oi_beams[0]->nr-1;
	print_info("Source plane %d is outside range\n",ir+1);
      }
      memset(nparr,0,par->kmap->num_pix*sizeof(int));

      for(ib=0;ib<par->n_beams_here;ib++) {
	int ind_cth;
	int icth0=par->oi_beams[ib]->icth0_arr[irb];
	int iphi0=par->oi_beams[ib]->iphi0_arr[irb];
	int ncth=par->oi_beams[ib]->icthf_arr[irb]-par->oi_beams[ib]->icth0_arr[irb]+1;
	int nphi=par->oi_beams[ib]->iphif_arr[irb]-par->oi_beams[ib]->iphi0_arr[irb]+1;
	int nside=par->oi_beams[ib]->nside_arr[irb];
	double dphi=M_PI/nside;
#if PIXTYPE==PT_CEA
	double dcth=2./nside;
#elif PIXTYPE==PT_CAR
	double dth=M_PI/nside;
#endif //PIXTYPE
	for(ind_cth=0;ind_cth<ncth;ind_cth++) {
	  int ind_phi;
	  int ind_cth_t=ind_cth*nphi;
#if PIXTYPE==PT_CEA
	  double cth0=get_cosine(ind_cth+icth0+0.0,dcth);
#elif PIXTYPE==PT_CAR
	  double cth0=get_cosine(ind_cth+icth0+0.0,dth);
	  double dcth=get_cosine(ind_cth+icth0+1.0,dth)-cth0;
#endif //PIXTYPE
	  int nsub_perside=(int)(sqrt(dcth*dphi*inv_hpix_area)+1);
	  double dcth_sub=dcth/nsub_perside;
	  double dphi_sub=dphi/nsub_perside;
	  for(ind_phi=0;ind_phi<nphi;ind_phi++) {
	    int icth_sub;
	    double phi0=(ind_phi+iphi0)*dphi;
	    double pxx=par->p_xx_beams[ib][irb][ind_cth_t+ind_phi];
	    double pyy=par->p_yy_beams[ib][irb][ind_cth_t+ind_phi];
	    double kappa=pxx+pyy;
	    for(icth_sub=0;icth_sub<nsub_perside;icth_sub++) {
	      int iphi_sub;
	      double cth=cth0+(icth_sub+0.5)*dcth_sub;
	      for(iphi_sub=0;iphi_sub<nsub_perside;iphi_sub++) {
		double phi=phi0+(iphi_sub+0.5)*dphi_sub;
		long ipix=he_ang2pix(par->kmap->nside,cth,phi);
		long pix_id=par->kmap->listpix[ipix];
		if(pix_id<0)
		  report_error(1,"NOOO\n");
		par->kmap->data[irad_t+pix_id]+=kappa;
		nparr[pix_id]++;
	      }
	    }
	  }
	}
      }

      for(ipx=0;ipx<par->kmap->num_pix;ipx++) {
      	if(nparr[ipx]>0)
      	  par->kmap->data[irad_t+ipx]/=nparr[ipx];
      }
    } //end omp for
    free(nparr);
  } //end omp parallel

  if(NodeThis==0) timer(2);
  printf("\n");
}
