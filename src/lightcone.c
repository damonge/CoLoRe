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

static inline void cart2sph(double x,double y,double z,double *r,double *cth,double *phi)
{
  *r=sqrt(x*x+y*y+z*z);

  if((*r)==0) {
    *cth=1;
    *phi=0;
  }
  else {
    double xn=x/(*r);
    double yn=y/(*r);

    *cth=z/(*r);
    if((xn==0)&&(yn==0))
      *phi=0;
    else {
      *phi=atan2(yn,xn);
      if((*phi)<0)
      	(*phi)+=2*M_PI;
    }
  }
}

static double get_rvel(ParamCoLoRe *par,int ix,int iy,int iz,
		       double x0,double y0,double z0,double rr)
{
  double v[3],u[3];
  double idx=par->n_grid/par->l_box;
  long ngx=2*(par->n_grid/2+1);
  long iz_hi=iz+1,iz_lo=iz-1,iz_0=iz;
  long iy_hi=iy+1,iy_lo=iy-1,iy_0=iy;
  long ix_hi=ix+1,ix_lo=ix-1,ix_0=ix;
  if(iy==0) iy_lo=par->n_grid-1;
  if(iy==par->n_grid-1) iy_hi=0;
  if(ix==0) ix_lo=par->n_grid-1;
  if(ix==par->n_grid-1) ix_hi=0;
  iz_0*=ngx*par->n_grid;
  iz_lo*=ngx*par->n_grid;
  iz_hi*=ngx*par->n_grid;
  iy_0*=ngx;
  iy_lo*=ngx;
  iy_hi*=ngx;

  u[0]=x0/rr; u[1]=y0/rr; u[2]=z0/rr;
  v[0]=par->grid_npot[ix_hi+iy_0+iz_0]-par->grid_npot[ix_lo+iy_0+iz_0];
  v[1]=par->grid_npot[ix_0+iy_hi+iz_0]-par->grid_npot[ix_0+iy_lo+iz_0];
  if(iz==0)
    v[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->slice_left[ix_0+iy_0];
  else if(iz==par->nz_here-1)
    v[2]=par->slice_right[ix_0+iy_0]-par->grid_npot[ix_0+iy_0+iz_lo];
  else
    v[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->grid_npot[ix_0+iy_0+iz_lo];

  return 0.5*idx*(v[0]*u[0]+v[1]*u[1]+v[2]*u[2]);
}

static void get_sources_cartesian_single(ParamCoLoRe *par,int ipop)
{
  int ii,nthr;
  int ngx=2*(par->n_grid/2+1);
  long *np_tot_thr;
  int *nsources=my_calloc(ngx*((long)(par->n_grid*par->nz_here)),sizeof(int));
#ifdef _HAVE_OMP
  nthr=omp_get_max_threads();
#else //_HAVE_OMP
  nthr=1;
#endif //_HAVE_OMP
  np_tot_thr=my_calloc(nthr,sizeof(long));

  print_info(" %d-th galaxy population\n",ipop);
  if(NodeThis==0) timer(0);
  print_info("   Poisson-sampling\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none)			\
  shared(par,np_tot_thr,IThread0,ipop,nthr,nsources)
#endif //_HAVE_OMP
  {
    long iz;
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    double dx=par->l_box/par->n_grid;
    double cell_vol=dx*dx*dx;
    int ngx=2*(par->n_grid/2+1);
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_srcs*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      double z0=(iz+par->iz0_here+0.5)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	double y0=(iy+0.5)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  int npp=0;
	  long index=ix+indexy+indexz;
	  double x0=(ix+0.5)*dx-par->pos_obs[0];
	  double r=sqrt(x0*x0+y0*y0+z0*z0);
	  double ndens=get_bg(par,r,BG_NZ_SRCS,ipop);
	  if(ndens>0) {
	    double bias=get_bg(par,r,BG_BZ_SRCS,ipop);
	    double dnorm=get_bg(par,r,BG_NORM_SRCS,ipop);
	    double lambda=ndens*cell_vol*bias_model(par->grid_dens[index],bias)*dnorm;
	    npp=rng_poisson(lambda,rng_thr);
	  }

	  nsources[index]=npp;
	  np_tot_thr[ithr]+=npp;
	}
      }
    }//end omp for

    end_rng(rng_thr);
  }//end omp parallel
  if(NodeThis==0) timer(2);

  par->nsources_this[ipop]=0;
  for(ii=0;ii<nthr;ii++)
    par->nsources_this[ipop]+=np_tot_thr[ii];

  long nsources_total=0;
#ifdef _HAVE_MPI
  MPI_Allreduce(&(par->nsources_this[ipop]),&nsources_total,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
#else //_HAVE_MPI
  nsources_total=par->nsources_this[ipop];
#endif //_HAVE_MPI

  print_info("   There will be %ld objects in total \n",(long)nsources_total);
#ifdef _DEBUG
  fprintf(par->f_dbg,"Node %d has %ld particles\n",NodeThis,(long)(par->nsources_this));
#endif //_DEBUG

  for(ii=nthr-1;ii>0;ii--) {
    int jj;
    long nh=0;
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
  shared(par,IThread0,np_tot_thr,ipop,nthr,nsources)
#endif //_HAVE_OMP
  {
    long iz;
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    double dx=par->l_box/par->n_grid;
    int ngx=2*(par->n_grid/2+1);
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_srcs*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);
    double factor_vel=-par->fgrowth_0/(1.5*par->hubble_0*par->OmegaM);

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      double z0=(iz+par->iz0_here+0.5)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	double y0=(iy+0.5)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  double x0=(ix+0.5)*dx-par->pos_obs[0];
	  long index=ix+indexy+indexz;
	  int npp=nsources[index];
	  if(npp>0) {
	    int ip;
	    double rr=sqrt(x0*x0+y0*y0+z0*z0);
	    double rvel=factor_vel*get_rvel(par,ix,iy,iz,x0,y0,z0,rr);
	    double dz_rsd=rvel*get_bg(par,rr,BG_V1,0);
	    for(ip=0;ip<npp;ip++) {
	      double cth,phi,r;
	      long pid=np_tot_thr[ithr];
	      double x=x0+dx*(rng_01(rng_thr)-0.5);
	      double y=y0+dx*(rng_01(rng_thr)-0.5);
	      double z=z0+dx*(rng_01(rng_thr)-0.5);
	      cart2sph(x,y,z,&r,&cth,&phi);
	      par->srcs[ipop][pid].ra=RTOD*phi;
	      par->srcs[ipop][pid].dec=90-RTOD*acos(cth);
	      par->srcs[ipop][pid].z0=get_bg(par,r,BG_Z,0);
	      par->srcs[ipop][pid].dz_rsd=dz_rsd;
	      par->srcs[ipop][pid].e1=-1;
	      par->srcs[ipop][pid].e2=-1;
	      np_tot_thr[ithr]++;
	    }
	  }
	}
      }
    }//end omp for
    end_rng(rng_thr);
  }//end omp parallel
  if(NodeThis==0) timer(2);

  free(np_tot_thr);
  free(nsources);
}  

static void get_sources_single(ParamCoLoRe *par,int ipop)
{
  //////
  // Uses the gaussian matter density field to obtain a
  // poisson sampling of point sources (returns an integer array
  // with the number of sources in each cell).

  int ii,nthr;
  long *np_tot_thr;
#ifdef _HAVE_OMP
  nthr=omp_get_max_threads();
#else //_HAVE_OMP
  nthr=1;
#endif //_HAVE_OMP
  np_tot_thr=my_calloc(nthr,sizeof(long));

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
      double ndens=get_bg(par,rm,BG_NZ_SRCS,ipop);
      if(ndens>0) {
	int ib;
	double bias=get_bg(par,rm,BG_BZ_SRCS,ipop);
	double dnorm=get_bg(par,rm,BG_NORM_SRCS,ipop);
	double pixarea=he_pixel_area(par->oi_beams[0]->nside_arr[ir]);
	double cell_vol=(rf*rf*rf-r0*r0*r0)*pixarea/3;
	for(ib=0;ib<par->n_beams_here;ib++) {
	  int ipix;
	  OnionInfo *oi=par->oi_beams[ib];
	  flouble *dens_slice=par->dens_beams[ib][ir];
	  int *nsrc_slice=par->nsrc_beams[ib][ir];
	  for(ipix=0;ipix<oi->num_pix[ir];ipix++) {
	    double lambda=ndens*cell_vol*bias_model(dens_slice[ipix],bias)*dnorm;
	    int npp=rng_poisson(lambda,rng_thr);
	    nsrc_slice[ipix]=npp;
	    np_tot_thr[ithr]+=npp;
	  }
	}
      }
    } //end omp for
    end_rng(rng_thr);
  } //end omp parallel

  par->nsources_this[ipop]=0;
  for(ii=0;ii<nthr;ii++)
    par->nsources_this[ipop]+=np_tot_thr[ii];

  long nsources_total=0;
#ifdef _HAVE_MPI
  MPI_Allreduce(&(par->nsources_this[ipop]),&nsources_total,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
#else //_HAVE_MPI
  nsources_total=par->nsources_this[ipop];
#endif //_HAVE_MPI

  print_info("   There will be %ld objects in total \n",(long)nsources_total);
#ifdef _DEBUG
  fprintf(par->f_dbg,"Node %d has %ld particles\n",NodeThis,(long)(par->nsources_this));
#endif //_DEBUG

  for(ii=nthr-1;ii>0;ii--) {
    int jj;
    long nh=0;
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
      double ndens=get_bg(par,rm,BG_NZ_SRCS,ipop);
      if(ndens>0) {
	int ib;
	for(ib=0;ib<par->n_beams_here;ib++) {
	  int ipix;
	  OnionInfo *oi=par->oi_beams[ib];
	  flouble *vrad_slice=par->vrad_beams[ib][ir];
	  int *nsrc_slice=par->nsrc_beams[ib][ir];
	  for(ipix=0;ipix<oi->num_pix[ir];ipix++) {
	    int ip;
	    double e1=0,e2=0;
	    int npp=nsrc_slice[ipix];
	    double dz_rsd=vrad_slice[ipix];
	    if(par->shear_srcs[ipop]) {
	      double pxx=par->p_xx_beams[ib][ir][ipix];
	      double pxy=par->p_xy_beams[ib][ir][ipix];
	      double pyy=par->p_yy_beams[ib][ir][ipix];
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
	      double th,phi;
	      double r=r0+dr*rng_01(rng_thr); //TODO: this is not completely correct
	      long pid=np_tot_thr[ithr];
	      get_random_angles(rng_thr,ipix,oi->ipix0_arr[ir],oi->nside_arr[ir],&th,&phi);
	      par->srcs[ipop][pid].ra=RTOD*phi;
	      par->srcs[ipop][pid].dec=90-RTOD*th;
	      par->srcs[ipop][pid].z0=get_bg(par,r,BG_Z,0);
	      par->srcs[ipop][pid].dz_rsd=dz_rsd;
	      par->srcs[ipop][pid].e1=e1;
	      par->srcs[ipop][pid].e2=e2;
	      np_tot_thr[ithr]++;
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
void integrate_isw(ParamCoLoRe *par)
{
  int ib;

#ifdef _DEBUG
  print_info("*** Integrating ISW\n");
  if(NodeThis==0) timer(0);
#endif //_DEBUG
  for(ib=0;ib<par->n_beams_here;ib++) {
    int ir;
    OnionInfo *oi=par->oi_beams[ib];
    int nside_old=oi->nside_arr[0];
    double *pdot_old=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pdot_new=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    for(ir=0;ir<oi->nr;ir++) {
      int ipix;
      int nside_new=oi->nside_arr[ir];
      int ratio_oldnew=nside_new/nside_old;
      double dr=oi->rf_arr[ir]-oi->r0_arr[ir];
      double g_phi=2*dr;
      for(ipix=0;ipix<oi->num_pix[ir];ipix++) {
	int ipix_old=ipix/(ratio_oldnew*ratio_oldnew);
	pdot_new[ipix]=pdot_old[ipix_old]+g_phi*par->pdot_beams[ib][ir][ipix];
	par->pdot_beams[ib][ir][ipix]=pdot_new[ipix];
      }
      nside_old=nside_new;
      memcpy(pdot_old,pdot_new,oi->num_pix[ir]*sizeof(double));
    }
    free(pdot_new); free(pdot_old);
  }
#ifdef _DEBUG
  if(NodeThis==0) timer(2);
#endif //_DEBUG
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
    OnionInfo *oi=par->oi_beams[ib];
    int nside_old=oi->nside_arr[0];
    double *pxx1_old=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pxx1_new=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pxx2_old=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pxx2_new=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pxy1_old=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pxy1_new=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pxy2_old=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pxy2_new=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pyy1_old=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pyy1_new=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pyy2_old=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    double *pyy2_new=my_calloc(oi->num_pix[oi->nr-1],sizeof(double));
    for(ir=0;ir<oi->nr;ir++) {
      int ipix;
      int nside_new=oi->nside_arr[ir];
      int ratio_oldnew=nside_new/nside_old;
      double r0=oi->r0_arr[ir];
      double rf=oi->rf_arr[ir];
      double rm=0.5*(r0+rf),dr=rf-r0;
      double integ1=rm*dr;
      double integ2=rm*rm*dr;
      //      double integ1=0.5*(rf*rf-r0*r0);
      //      double integ2=(rf*rf*rf-r0*r0*r0)/3;
      for(ipix=0;ipix<oi->num_pix[ir];ipix++) {
	int ipix_old=ipix/(ratio_oldnew*ratio_oldnew);
	pxx1_new[ipix]=pxx1_old[ipix_old]+integ1*par->p_xx_beams[ib][ir][ipix];
	pxx2_new[ipix]=pxx2_old[ipix_old]+integ2*par->p_xx_beams[ib][ir][ipix];
	par->p_xx_beams[ib][ir][ipix]=pxx1_new[ipix]-pxx2_new[ipix]/rf;
	pxy1_new[ipix]=pxy1_old[ipix_old]+integ1*par->p_xy_beams[ib][ir][ipix];
	pxy2_new[ipix]=pxy2_old[ipix_old]+integ2*par->p_xy_beams[ib][ir][ipix];
	par->p_xy_beams[ib][ir][ipix]=pxy1_new[ipix]-pxy2_new[ipix]/rf;
	pyy1_new[ipix]=pyy1_old[ipix_old]+integ1*par->p_yy_beams[ib][ir][ipix];
	pyy2_new[ipix]=pyy2_old[ipix_old]+integ2*par->p_yy_beams[ib][ir][ipix];
	par->p_yy_beams[ib][ir][ipix]=pyy1_new[ipix]-pyy2_new[ipix]/rf;
      }
      nside_old=nside_new;
      memcpy(pxx1_old,pxx1_new,oi->num_pix[ir]*sizeof(double));
      memcpy(pxx2_old,pxx2_new,oi->num_pix[ir]*sizeof(double));
      memcpy(pxy1_old,pxy1_new,oi->num_pix[ir]*sizeof(double));
      memcpy(pxy2_old,pxy2_new,oi->num_pix[ir]*sizeof(double));
      memcpy(pyy1_old,pyy1_new,oi->num_pix[ir]*sizeof(double));
      memcpy(pyy2_old,pyy2_new,oi->num_pix[ir]*sizeof(double));
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
  for(ipop=0;ipop<par->n_srcs;ipop++) {
    if(par->need_onions)
      get_sources_single(par,ipop);
    else
      get_sources_cartesian_single(par,ipop);
  }
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
  shell->num_pix=0;
  shell->listpix=my_malloc(npx*sizeof(long));

  for(ip=0;ip<npx;ip++) {
    int ib;
    int goodpix=0;

    if(par->need_onions) {
      int nside_ratio=shell->nside/par->oi_beams[0]->nside_arr[0];
      long ip_base=ip/(nside_ratio*nside_ratio);
      for(ib=0;ib<par->n_beams_here;ib++) {
	if(ip_base==par->oi_beams[ib]->ipix0_arr[0])
	  goodpix=1;
      }
    }
    else
      goodpix=1;
    
    if(goodpix) {
      shell->listpix[ip]=shell->num_pix;
      shell->num_pix++;
    }
    else
      shell->listpix[ip]=-1;
  }


  shell->data=my_calloc(shell->nr*shell->num_pix,sizeof(flouble));
  shell->nadd=my_calloc(shell->nr*shell->num_pix,sizeof(int));
}

static int *get_nsub_r_list(OnionInfo *oi,HealpixShells *sh)
{
  int ir;
  int irad=0;
  int *nr_sub=my_malloc(oi->nr*sizeof(int));

  for(ir=0;ir<oi->nr;ir++) {
    double rm=0.5*(oi->rf_arr[ir]+oi->r0_arr[ir]);
    irad=get_r_index_imap(sh,rm,irad);
    if((irad>=0) && (irad<sh->nr)) {
      double dr_hp=sh->rf[irad]-sh->r0[irad];
      double dr_pix=(oi->rf_arr[ir]-oi->r0_arr[ir]);
      nr_sub[ir]=(int)(MAX(dr_hp/dr_pix,1.));
    }
    else
      nr_sub[ir]=1;
  }

  return nr_sub;
}

static double get_rmax(HealpixShells *sh)
{
  int ir;
  double rmax=sh->rf[0];

  for(ir=0;ir<sh->nr;ir++) {
    if(sh->rf[ir]>rmax)
      rmax=sh->rf[ir];
  }

  return rmax;
}

static double get_min_vol(HealpixShells *sh)
{
  int ir;
  double pixel_area=he_pixel_area(sh->nside);
  double volmin=1E60;

  for(ir=0;ir<sh->nr;ir++) {
    double r0=sh->r0[ir];
    double rf=sh->rf[ir];
    double vol=pixel_area*(rf*rf*rf-r0*r0*r0);
    if(vol<volmin)
      volmin=vol;
  }

  return volmin;
}

static void get_imap_cartesian_single(ParamCoLoRe *par,int ipop)
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
    long iz;
    double dx=par->l_box/par->n_grid;
    int ngx=2*(par->n_grid/2+1);
    double factor_vel=-par->fgrowth_0/(1.5*par->hubble_0*par->OmegaM);
    double hpa_third=4*M_PI/(3.*he_nside2npix(par->imap[ipop]->nside));

    double volmin=get_min_vol(par->imap[ipop]);
    int n_resample=(int)(MAX((8*dx*dx*dx/volmin),1.));
    double *disp=my_calloc(n_resample*3,sizeof(double));
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_imap*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);
    for(iz=0;iz<3*n_resample;iz++)
      disp[iz]=rng_01(rng_thr);
    end_rng(rng_thr);

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      int irad=0;
      long indexz=iz*((long)(ngx*par->n_grid));
      double z0=(iz+par->iz0_here+0.5)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	double y0=(iy+0.5)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  long index=ix+indexy+indexz;
	  double x0=(ix+0.5)*dx-par->pos_obs[0];
	  double r0=sqrt(x0*x0+y0*y0+z0*z0);
	  double tmean=get_bg(par,r0,BG_TZ_IMAP,ipop);
	  if(tmean>0) {
	    int isub,n_resample_here;
	    double bias=get_bg(par,r0,BG_BZ_IMAP,ipop);
	    double dnorm=get_bg(par,r0,BG_NORM_IMAP,ipop);
	    double rvel=factor_vel*get_rvel(par,ix,iy,iz,x0,y0,z0,r0);
	    double dr_rsd=rvel*get_bg(par,r0,BG_V1,0)*get_bg(par,r0,BG_IH,0);
	    double temp=tmean*bias_model(par->grid_dens[index],bias)*dnorm;
	    if((irad>=0) && (irad<par->imap[ipop]->nr)) {
	      double pixvol;
	      double rfh,r0h;
	      irad=get_r_index_imap(par->imap[ipop],r0,irad);
	      rfh=par->imap[ipop]->rf[irad];
	      r0h=par->imap[ipop]->r0[irad];
	      pixvol=(rfh*rfh*rfh-r0h*r0h*r0h)*hpa_third;
	      n_resample_here=(int)(MAX((8*dx*dx*dx/pixvol),1.));
	    }
	    else
	      n_resample_here=n_resample;
	    if(n_resample_here>n_resample) {
	      printf("SHIT\n");
	      exit(1);
	    }
	    n_resample_here=n_resample;
	    for(isub=0;isub<n_resample_here;isub++) {
	      double cth,phi,r;
	      double x=x0+(disp[3*isub+0]-0.5)*dx;
	      double y=y0+(disp[3*isub+1]-0.5)*dx;
	      double z=z0+(disp[3*isub+2]-0.5)*dx;
	      cart2sph(x,y,z,&r,&cth,&phi);
	      irad=get_r_index_imap(par->imap[ipop],r+dr_rsd,irad);
	      if((irad>=0) && (irad<par->imap[ipop]->nr)) {
		long pix_id;
		long irad_t=irad*par->imap[ipop]->num_pix;
		long pix_id_ring=he_ang2pix(par->imap[ipop]->nside,cth,phi);
		ring2nest(par->imap[ipop]->nside,pix_id_ring,&pix_id);
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
		par->imap[ipop]->data[irad_t+pix_id]+=temp;
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
		par->imap[ipop]->nadd[irad_t+pix_id]++;
	      }
	    }
	  }
	}
      }
    }//end omp for
    free(disp);

  }//end omp parallel
  if(NodeThis==0) timer(2);
}

static void get_imap_single(ParamCoLoRe *par,int ipop)
{
  print_info(" %d-th IM species\n",ipop);
  if(NodeThis==0) timer(0);
  find_shell_pixels(par,par->imap[ipop]);

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,IThread0,ipop)
#endif //_HAVE_OMP
  {
    int ir;
    double rmax_imap=get_rmax(par->imap[ipop])+20.; //We add 20 Mpc/h extra for RSDs
    int *nsub_r_arr=get_nsub_r_list(par->oi_beams[0],par->imap[ipop]);

#ifdef _HAVE_OMP
#pragma omp for schedule(dynamic)
#endif //_HAVE_OMP
    for(ir=0;ir<par->oi_beams[0]->nr;ir++) {
      double rf=par->oi_beams[0]->rf_arr[ir];
      if(rf<rmax_imap) {
	double r0=par->oi_beams[0]->r0_arr[ir];
	double rm=(rf+r0)*0.5;
	double dr=rf-r0;
	double dr_sub=dr/nsub_r_arr[ir];
	double tmean=get_bg(par,rm,BG_TZ_IMAP,ipop);
	if(tmean>0) {
	  int ib;
	  int irad=0;
	  double bias=get_bg(par,rm,BG_BZ_IMAP,ipop);
	  double dnorm=get_bg(par,rm,BG_NORM_IMAP,ipop);
	  double prefac_rsd=get_bg(par,rm,BG_IH,0);
	  int nsub_perside=MAX((par->imap[ipop]->nside/par->oi_beams[0]->nside_arr[ir]),1);
	  int nside_ratio=MAX((par->oi_beams[0]->nside_arr[ir]/par->imap[ipop]->nside),1);
	  for(ib=0;ib<par->n_beams_here;ib++) {
	    int ipix;
	    OnionInfo *oi=par->oi_beams[ib];
	    flouble *dens_slice=par->dens_beams[ib][ir];
	    flouble *vrad_slice=par->vrad_beams[ib][ir];
	    for(ipix=0;ipix<oi->num_pix[ir];ipix++) {
	      int ipix_sub;
	      double temp=tmean*bias_model(dens_slice[ipix],bias)*dnorm;
	      double dr_rsd=prefac_rsd*vrad_slice[ipix];
	      long ipix0=(oi->ipix0_arr[ir]+ipix)*nsub_perside*nsub_perside;
	      for(ipix_sub=0;ipix_sub<nsub_perside*nsub_perside;ipix_sub++) {
		int ir_sub;
		long pix_id,ip_shell=(ipix0+ipix_sub)/(nside_ratio*nside_ratio);
		pix_id=par->imap[ipop]->listpix[ip_shell];
		if(pix_id<0)
		  report_error(1,"NOOOO \n");
		for(ir_sub=0;ir_sub<nsub_r_arr[ir];ir_sub++) {
		  double r=r0+dr_rsd+dr_sub*ir_sub;
		  irad=get_r_index_imap(par->imap[ipop],r,irad);
		  if((irad>=0) && (irad<par->imap[ipop]->nr)) {
		    long irad_t=irad*par->imap[ipop]->num_pix;
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
		    par->imap[ipop]->data[irad_t+pix_id]+=temp;
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
                    par->imap[ipop]->nadd[irad_t+pix_id]++;
                  }
		}
	      }
	    }
	  }
	}
      }
    } //end omp for
    free(nsub_r_arr);
  } //end omp parallel

  if(NodeThis==0) timer(2);
}

void get_imap(ParamCoLoRe *par)
{
  int ipop;

  //First, compute lensing Hessian
  print_info("*** Getting intensity maps\n");
  for(ipop=0;ipop<par->n_imap;ipop++) {
    if(par->need_onions)
      get_imap_single(par,ipop);
    else
      get_imap_cartesian_single(par,ipop);
  }
  print_info("\n");
}

void get_kappa(ParamCoLoRe *par)
{

  print_info("*** Getting kappa maps\n");
  if(NodeThis==0) timer(0);
  find_shell_pixels(par,par->kmap);

#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par)
#endif //_HAVE_OMP
  {
    int ir;

    //Maybe OMP this
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ir=0;ir<par->n_kappa;ir++) {
      int ib,nsub_perside,nside_ratio,irb=0;
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

      nsub_perside=MAX((par->kmap->nside/par->oi_beams[0]->nside_arr[irb]),1);
      nside_ratio=MAX((par->oi_beams[0]->nside_arr[irb]/par->kmap->nside),1);
      for(ib=0;ib<par->n_beams_here;ib++) {
	int ipix;
	OnionInfo *oi=par->oi_beams[ib];
	for(ipix=0;ipix<oi->num_pix[irb];ipix++) {
	  int ipix_sub;
	  double pxx=par->p_xx_beams[ib][irb][ipix];
	  double pyy=par->p_yy_beams[ib][irb][ipix];
	  double kappa=pxx+pyy;
	  long ipix0=(oi->ipix0_arr[irb]+ipix)*nsub_perside*nsub_perside;
	  for(ipix_sub=0;ipix_sub<nsub_perside*nsub_perside;ipix_sub++) {
	    long pix_id,ip_shell=(ipix0+ipix_sub)/(nside_ratio*nside_ratio);
	    pix_id=par->kmap->listpix[ip_shell];
	    if(pix_id<0)
	      report_error(1,"NOOOO \n");
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
	    par->kmap->data[irad_t+pix_id]+=kappa;
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
	    par->kmap->nadd[irad_t+pix_id]++;
	  }
	}
      }
    } //end omp for
  } //end omp parallel

  if(NodeThis==0) timer(2);
  print_info("\n");
}

void get_isw(ParamCoLoRe *par)
{

  print_info("*** Getting isw maps\n");
  if(NodeThis==0) timer(0);
  find_shell_pixels(par,par->pd_map);

#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par)
#endif //_HAVE_OMP
  {
    int ir;

    //Maybe OMP this
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ir=0;ir<par->n_isw;ir++) {
      int ib,nsub_perside,nside_ratio,irb=0;
      double r=par->pd_map->r0[ir];
      long irad_t=ir*par->pd_map->num_pix;
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

      nsub_perside=MAX((par->pd_map->nside/par->oi_beams[0]->nside_arr[irb]),1);
      nside_ratio=MAX(((par->oi_beams[0]->nside_arr[irb]*nsub_perside)/par->pd_map->nside),1);
      for(ib=0;ib<par->n_beams_here;ib++) {
	long ipix;
	OnionInfo *oi=par->oi_beams[ib];
	for(ipix=0;ipix<oi->num_pix[irb];ipix++) {
	  int ipix_sub;
	  double isw=par->pdot_beams[ib][irb][ipix];
	  long ipix0=(oi->ipix0_arr[irb]+ipix)*nsub_perside*nsub_perside;
	  for(ipix_sub=0;ipix_sub<nsub_perside*nsub_perside;ipix_sub++) {
	    long pix_id,ip_shell=(ipix0+ipix_sub)/(nside_ratio*nside_ratio);
	    pix_id=par->pd_map->listpix[ip_shell];
	    if(pix_id<0)
	      report_error(1,"NOOOO \n");
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
	    par->pd_map->data[irad_t+pix_id]+=isw;
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
	    par->pd_map->nadd[irad_t+pix_id]++;
	  }
	}
      }
    } //end omp for
  } //end omp parallel

  if(NodeThis==0) timer(2);
  print_info("\n");
}
