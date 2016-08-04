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

inline void cart2sph(double x,double y,double z,double *r,double *cth,double *phi)
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

static void get_sources_single(ParamCoLoRe *par,int ipop,int *nsources)
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
#pragma omp parallel default(none)			\
  shared(par,np_tot_thr,IThread0,ipop,nthr,nsources)
#endif //_HAVE_OMP
  {
    lint iz;
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    double dx=par->l_box/par->n_grid;
    double cell_vol=dx*dx*dx;
    int ngx=2*(par->n_grid/2+1);
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_gals*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint indexz=iz*((lint)(ngx*par->n_grid));
      double z0=(iz+par->iz0_here+0.5)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint indexy=iy*ngx;
	double y0=(iy+0.5)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  int npp=0;
	  lint index=ix+indexy+indexz;
	  double x0=(ix+0.5)*dx-par->pos_obs[1];
	  double r=sqrt(x0*x0+y0*y0+z0*z0);
	  double redshift=z_of_r(par,r);
	  double ndens=ndens_of_z_gals(par,redshift,ipop);
	  if(ndens>0) {
	    double bias=bias_of_z_gals(par,redshift,ipop);
	    double gfb=dgrowth_of_r(par,r)*bias;
	    double lambda=ndens*cell_vol*
	      exp(gfb*(par->grid_dens[index]-0.5*gfb*par->sigma2_gauss));
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

  lint nsources_total=0;
#ifdef _HAVE_MPI
  MPI_Allreduce(&(par->nsources_this[ipop]),&nsources_total,1,LINT_MPI,MPI_SUM,MPI_COMM_WORLD);
#else //_HAVE_MPI
  nsources_total=par->nsources_this[ipop];
#endif //_HAVE_MPI

  print_info("   There will be %ld objects in total \n",(long)nsources_total);
  //#ifdef _DEBUG
  //  printf("Node %d has %ld particles\n",NodeThis,(long)(par->nsources_this));
  //#endif //_DEBUG

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
  shared(par,IThread0,np_tot_thr,ipop,nthr,nsources)
#endif //_HAVE_OMP
  {
    lint iz;
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    double dx=par->l_box/par->n_grid;
    int ngx=2*(par->n_grid/2+1);
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_gals*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint indexz=iz*((lint)(ngx*par->n_grid));
      double z0=(iz+par->iz0_here)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint indexy=iy*ngx;
	double y0=iy*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  double x0=ix*dx-par->pos_obs[1];
	  lint index=ix+indexy+indexz;
	  int npp=nsources[index];
	  if(npp>0) {
	    int ip;
	    double rr=sqrt(x0*x0+y0*y0+z0*z0);
	    double dz_rsd=par->grid_rvel[index]*vgrowth_of_r(par,rr);
	    for(ip=0;ip<npp;ip++) {
	      double cth,phi,r;
	      lint pid=np_tot_thr[ithr];
	      double x=x0+dx*rng_01(rng_thr);
	      double y=y0+dx*rng_01(rng_thr);
	      double z=z0+dx*rng_01(rng_thr);
	      cart2sph(x,y,z,&r,&cth,&phi);
	      par->gals[ipop][pid].ra=RTOD*phi;
	      par->gals[ipop][pid].dec=90-RTOD*acos(cth);
	      par->gals[ipop][pid].z0=z_of_r(par,r);
	      par->gals[ipop][pid].dz_rsd=dz_rsd;
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
}

void get_psi_potential(ParamCoLoRe *par)
{
  int ix, iy, iz;
  int iplane;
  int iplane1;
  long ipix;
  int ngx=2*(par->n_grid/2+1);
  int npix = 12*(par->nside)*(par->nside);
  int nplanes = par->n_lens_planes;
  double dx=par->l_box/par->n_grid;
  int imem;
  flouble **phi_potential = (flouble**)malloc(nplanes*sizeof(flouble*));
  for (imem=0; imem<nplanes; imem++){
    phi_potential[imem]=my_malloc(npix*sizeof(flouble));
  }
  if(NodeThis==0) timer(0);
  print_info("   Computing lensing potential\n");
    #ifdef _HAVE_OMP
  #pragma omp parallel default(none)			\
    shared(par, ix, iy, iz, iplane, iplane1, ipix, phi_potential, dx)
  #endif //_HAVE_OMP
    {
  #ifdef _HAVE_OMP
  #pragma omp for schedule(static) collapse(3)
  #endif //_HAVE_OMP
  for(ix=0; ix<par->n_grid; ix++){
    for(iy=0; iy<par->n_grid; iy++){
      for(iz=0; iz<par->n_grid; iz++){
        iplane = (int) sqrt(ix*ix+iy*iy+iz*iz)/sqrt(3.*par->n_grid*par->n_grid)*nplanes;
        double myvec[3];
        myvec[0]=ix-par->n_grid/2;
        myvec[1]=iy-par->n_grid/2;
        myvec[2]=iz-par->n_grid/2;
        vec2pix_ring(par->nside,myvec,&ipix);
        double x0=(ix+0.5)*dx-par->pos_obs[1];
        double y0=iy*dx-par->pos_obs[1];
        double z0=(iz+par->iz0_here)*dx-par->pos_obs[2];
        double r=sqrt(x0*x0+y0*y0+z0*z0);
        lint indexz=iz*((lint)(ngx*par->n_grid));
        lint indexy=iy*ngx;
        lint index=ix+indexy+indexz;
        phi_potential[iplane][ipix]=-3/2*par->hubble_0*par->OmegaM/par->fgrowth_0*dgrowth_of_r(par,r)*par->grid_vpot[index]*KMTOMPC;
      }
    }
  }
  par->psi_potential = my_calloc(nplanes,npix*sizeof(flouble));
  #ifdef _HAVE_OMP
  #pragma omp for schedule(static) collapse(3)
  #endif //_HAVE_OMP
  for(ipix=0; ipix<npix; ipix++){
    for(iplane=0; iplane<nplanes; iplane++){
      for(iplane1=0; iplane1<iplane; iplane1++){
        lint index = iplane*npix+ipix;
        //Only implemented for Flat FRWL
        par->psi_potential[index]+=-2*par->l_box*(iplane-iplane1)/nplanes*phi_potential[iplane1][ipix];
      }
    }
  }
}
  if(NodeThis==0) timer(2);
  free(phi_potential);
}


void get_sources(ParamCoLoRe *par)
{
  int ipop;
  int ngx=2*(par->n_grid/2+1);
  int *nsources=my_malloc(ngx*((lint)(par->n_grid*par->n_grid))*sizeof(int));

  print_info("*** Getting point sources\n");
  for(ipop=0;ipop<par->n_gals;ipop++)
    get_sources_single(par,ipop,nsources);

  free(nsources);
  print_info("\n");
}
