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
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_gals*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);

#ifdef _HAVE_OMP
#pragma omp for schedule(static) //TODO: this will give very bad load balance
#endif //_HAVE_OMP 
    for(ir=0;ir<par->oi_beams[0]->nr;ir++) {
      double r0=par->oi_beams[0]->r0_arr[ir];
      double rf=par->oi_beams[0]->rf_arr[ir];
      double rm=(rf+r0)*0.5;
      double redshift=z_of_r(par,rm);
      double ndens=ndens_of_z_gals(par,redshift,ipop);
      if(ndens>0) {
	int ib;
	int nside=par->oi_beams[0]->nside_arr[ir];
	double dcth=2./nside;
	double dphi=M_PI/nside;
	double cell_vol=(rf*rf*rf-r0*r0*r0)*dcth*dphi/3;
	double bias=bias_of_z_gals(par,redshift,ipop);
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
  
  par->gals[ipop]=my_malloc(par->nsources_this[ipop]*sizeof(Gal));

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
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_gals*IThread0);
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
      double ndens=ndens_of_z_gals(par,redshift,ipop);
      if(ndens>0) {
	int ib;
	int nside=par->oi_beams[0]->nside_arr[ir];
	double dcth=2./nside;
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
	    double cth0=(ind_cth+icth0)*dcth-1;
	    for(ind_phi=0;ind_phi<nphi;ind_phi++) {
	      int ip;
	      double e1=0,e2=0;
	      double phi0=(ind_phi+iphi0)*dphi;
	      int npp=nsrc_slice[ind_cth_t+ind_phi];
	      double dz_rsd=vg*vrad_slice[ind_cth_t+ind_phi];
	      if(par->do_lensing) {
		double pxx=par->p_xx_beams[ib][ir][ind_cth_t+ind_phi];
		double pxy=par->p_xy_beams[ib][ir][ind_cth_t+ind_phi];
		double pyy=par->p_yy_beams[ib][ir][ind_cth_t+ind_phi];
		double kappa=-(pxx+pyy)*0.5;
		double g1=-(pxx-pyy)*0.5;
		double g2=-pxy;
		double fac=2*(1-kappa)/((1-kappa)*(1-kappa)+g1*g1+g2*g2);
		e1=fac*g1;
		e2=fac*g2;
	      }
	      for(ip=0;ip<npp;ip++) {
		double cth=cth0+dcth*rng_01(rng_thr);
		double phi=phi0+dphi*rng_01(rng_thr);
		double r=r0+dr*rng_01(rng_thr); //TODO: this is not completely correct
		lint pid=np_tot_thr[ithr];
		par->gals[ipop][pid].ra=RTOD*phi;
		par->gals[ipop][pid].dec=90-RTOD*acos(cth);
		par->gals[ipop][pid].z0=z_of_r(par,r);
		par->gals[ipop][pid].dz_rsd=dz_rsd;
		par->gals[ipop][pid].e1=e1;
		par->gals[ipop][pid].e2=e2;
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

void integrate_lensing(ParamCoLoRe *par)
{
  int ib;

#ifdef _DEBUG
  print_info("Integrating lensing\n");
  if(NodeThis==0) timer(0);
#endif //_DEBUG
  for(ib=0;ib<par->n_beams_here;ib++) {
    int ir;
    double nside_old=par->oi_beams[ib]->nside_arr[0];
    double nphi_old=par->oi_beams[ib]->iphif_arr[0]-par->oi_beams[ib]->iphi0_arr[0];
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
      int ncth=par->oi_beams[ib]->icthf_arr[ir]-par->oi_beams[ib]->icth0_arr[ir];
      int nphi=par->oi_beams[ib]->iphif_arr[ir]-par->oi_beams[ib]->iphi0_arr[ir];
      double r0=par->oi_beams[ib]->r0_arr[ir];
      double rf=par->oi_beams[ib]->rf_arr[ir];
      double rm=0.5*(r0+rf);
      double redshift=z_of_r(par,rm);
      double g_phi=dgrowth_of_r(par,rm)*(1+redshift);
      int nside_new=par->oi_beams[ib]->nside_arr[ir];
      int nside_ratio=nside_new/nside_old;
      for(icth=0;icth<ncth;icth++) {
	int iphi;
	int icth_old=icth/nside_ratio;
	for(iphi=0;iphi<nphi;iphi++) {
	  int iphi_old=iphi/nside_ratio;
	  int index_new=iphi+nphi*icth;
	  int index_old=iphi_old+nphi_old*icth_old;
	  pxx1_new[index_new]=pxx1_old[index_old]+g_phi*par->p_xx_beams[ib][ir][index_new]*(rf*rf-r0*r0)/2;
	  pxx2_new[index_new]=pxx1_old[index_old]+g_phi*par->p_xx_beams[ib][ir][index_new]*(rf*rf*rf-r0*r0*r0)/3;
	  par->p_xx_beams[ib][ir][index_new]=2*(pxx1_new[index_new]-pxx1_new[index_new]/rf);
	  pxy1_new[index_new]=pxy1_old[index_old]+g_phi*par->p_xy_beams[ib][ir][index_new]*(rf*rf-r0*r0)/2;
	  pxy2_new[index_new]=pxy1_old[index_old]+g_phi*par->p_xy_beams[ib][ir][index_new]*(rf*rf*rf-r0*r0*r0)/3;
	  par->p_xy_beams[ib][ir][index_new]=2*(pxy1_new[index_new]-pxy1_new[index_new]/rf);
	  pyy1_new[index_new]=pyy1_old[index_old]+g_phi*par->p_yy_beams[ib][ir][index_new]*(rf*rf-r0*r0)/2;
	  pyy2_new[index_new]=pyy1_old[index_old]+g_phi*par->p_yy_beams[ib][ir][index_new]*(rf*rf*rf-r0*r0*r0)/3;
	  par->p_yy_beams[ib][ir][index_new]=2*(pyy1_new[index_new]-pyy1_new[index_new]/rf);
	}
      }
      nphi_old=nphi;
      nside_old=nside_new;
      memcpy(pxx1_old,pxx1_new,par->oi_beams[ib]->num_pix[ir]*sizeof(flouble));
      memcpy(pxx2_old,pxx2_new,par->oi_beams[ib]->num_pix[ir]*sizeof(flouble));
      memcpy(pxy1_old,pxy1_new,par->oi_beams[ib]->num_pix[ir]*sizeof(flouble));
      memcpy(pxy2_old,pxy2_new,par->oi_beams[ib]->num_pix[ir]*sizeof(flouble));
      memcpy(pyy1_old,pyy1_new,par->oi_beams[ib]->num_pix[ir]*sizeof(flouble));
      memcpy(pyy2_old,pyy2_new,par->oi_beams[ib]->num_pix[ir]*sizeof(flouble));
    }
    free(pxx1_old);
    free(pxx2_old);
    free(pxy1_old);
    free(pxy2_old);
    free(pyy1_old);
    free(pyy2_old);
  }
#ifdef _DEBUG
  if(NodeThis==0) timer(2);
#endif //_DEBUG
}

void get_sources(ParamCoLoRe *par)
{
  int ipop;

  if(par->do_lensing)
    integrate_lensing(par);

  print_info("*** Getting point sources\n");
  for(ipop=0;ipop<par->n_gals;ipop++)
    get_sources_single(par,ipop);
  print_info("\n");
}
