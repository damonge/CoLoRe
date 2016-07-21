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

static void compute_sigma_dens(ParamCoLoRe *par)
{
  double mean_gauss=0;
  par->sigma2_gauss=0;

  //Compute Gaussian variance
#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par,mean_gauss)
#endif //_HAVE_OMP
  {
    int iz;
    double sigma2_thr=0;
    double mean_thr=0;
    lint ng_tot=par->n_grid*((lint)(par->n_grid*par->n_grid));

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint iz0=iz*((lint)(2*(par->n_grid/2+1)*par->n_grid));
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint iy0=iy*2*(par->n_grid/2+1);
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=ix+iy0+iz0;
	  sigma2_thr+=par->grid_dens[index]*par->grid_dens[index];
	  mean_thr+=par->grid_dens[index];
	}
      }
    } //end omp for
#ifdef _HAVE_OMP
#pragma omp critical
#endif //_HAVE_OMP
    {
      mean_gauss+=mean_thr/ng_tot;
      par->sigma2_gauss+=sigma2_thr/ng_tot;
    } //end omp critical
  } //end omp parallel

#ifdef _HAVE_MPI
  double sigma2_gathered=0,mean_gathered=0;

  MPI_Allreduce(&(par->sigma2_gauss),&sigma2_gathered,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&mean_gauss,&mean_gathered,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  mean_gauss=mean_gathered;
  par->sigma2_gauss=sigma2_gathered;
#endif //_HAVE_MPI

  par->sigma2_gauss-=mean_gauss*mean_gauss;
  print_info(" <d>=%.3lE, <d^2>=%.3lE\n",mean_gauss,sqrt(par->sigma2_gauss));
}

static void fftw_wrap(int ng,dftw_complex *pin,flouble *pout)
{
#ifdef _SPREC
  fftwf_plan plan_ft;
#ifdef _HAVE_MPI
  plan_ft=fftwf_mpi_plan_dft_c2r_3d(ng,ng,ng,pin,pout,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else //_HAVE_MPI
  plan_ft=fftwf_plan_dft_c2r_3d(ng,ng,ng,pin,pout,FFTW_ESTIMATE);
#endif //_HAVE_MPI
  fftwf_execute(plan_ft);
  fftwf_destroy_plan(plan_ft);
#else //_SPREC
  fftw_plan plan_ft;
#ifdef _HAVE_MPI
  plan_ft=fftw_mpi_plan_dft_c2r_3d(ng,ng,ng,pin,pout,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else //_HAVE_MPI
  plan_ft=fftw_plan_dft_c2r_3d(ng,ng,ng,pin,pout,FFTW_ESTIMATE);
#endif //_HAVE_MPI
  fftw_execute(plan_ft);
  fftw_destroy_plan(plan_ft);
#endif //_SPREC
}

void init_fftw(ParamCoLoRe *par)
{
  ptrdiff_t dsize;
#ifdef _HAVE_MPI
  ptrdiff_t nz,iz0;

  //Initialize OpenMP fftw (if possible)
#ifdef _HAVE_OMP
  if(MPIThreadsOK) {
    int stat;
#ifdef _SPREC
    stat=fftwf_init_threads();
#else //_SPREC
    stat=fftw_init_threads();
#endif //_SPREC
    if(!stat) {
      fprintf(stderr,"Couldn't initialize FFTW threads \n");
      exit(1);
    }
  }
#endif //_HAVE_OMP

  //Initialize MPI fftw
#ifdef _SPREC
  fftwf_mpi_init();
#else //_SPREC
  fftw_mpi_init();
#endif //_SPREC

  //Plan OpenMP fftw (if possible)
#ifdef _HAVE_OMP
  if(MPIThreadsOK) {
#ifdef _SPREC
    fftwf_plan_with_nthreads(omp_get_max_threads());
#else //_SPREC
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif //_SPREC
  }
#endif //_HAVE_OMP

  //Get MPI FFT bounds
#ifdef _SPREC
  dsize=fftwf_mpi_local_size_3d(par->n_grid,par->n_grid,par->n_grid/2+1,MPI_COMM_WORLD,&nz,&iz0);
#else //_SPREC
  dsize=fftw_mpi_local_size_3d(par->n_grid,par->n_grid,par->n_grid/2+1,MPI_COMM_WORLD,&nz,&iz0);
#endif //_SPREC
  par->nz_here=nz;
  par->iz0_here=iz0;

#else //_HAVE_MPI
  int stat;
#ifdef _HAVE_OMP
#ifdef _SPREC
  stat=fftwf_init_threads();
#else //_SPREC
  stat=fftw_init_threads();
#endif //_SPREC
  if(!stat) {
    fprintf(stderr,"Couldn't initialize FFTW threads \n");
    exit(1);
  }
#ifdef _SPREC
  fftwf_plan_with_nthreads(omp_get_max_threads());
#else //_SPREC
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif //_SPREC
#endif //_HAVE_OMP

  dsize=(par->n_grid/2+1)*((lint)(par->n_grid*par->n_grid));
  par->nz_here=par->n_grid;
  par->iz0_here=0;
#endif //_HAVE_MPI

#ifdef _SPREC
  par->grid_dens_f=fftwf_alloc_complex(dsize);
#else //_SPREC
  par->grid_dens_f=fftw_alloc_complex(dsize);
#endif //_SPREC
  if(par->grid_dens_f==NULL)
    report_error(1,"Ran out of memory\n");
  par->grid_dens=(flouble *)(par->grid_dens_f);

#ifdef _SPREC
  par->grid_vpot_f=fftwf_alloc_complex(dsize);
#else //_SPREC
  par->grid_vpot_f=fftw_alloc_complex(dsize);
#endif //_SPREC
  if(par->grid_vpot_f==NULL)
    report_error(1,"Ran out of memory\n");
  par->grid_vpot=(flouble *)(par->grid_vpot_f);

  #ifdef _SPREC
    par->grid_npot_f=fftwf_alloc_complex(dsize);
  #else //_SPREC
    par->grid_npot_f=fftw_alloc_complex(dsize);
  #endif //_SPREC
    if(par->grid_npot_f==NULL)
      report_error(1,"Ran out of memory\n");
    par->grid_npot=(flouble *)(par->grid_npot_f);

  par->grid_rvel=my_malloc(2*dsize*sizeof(flouble));
#ifdef _HAVE_MPI
  par->slice_left=my_malloc(2*(par->n_grid/2+1)*par->n_grid*sizeof(flouble));
  par->slice_right=my_malloc(2*(par->n_grid/2+1)*par->n_grid*sizeof(flouble));
#endif //_HAVE_MPI
}

void end_fftw(void)
{
#ifdef _HAVE_MPI

#ifdef _HAVE_OMP
  if(MPIThreadsOK) {
#ifdef _SPREC
    fftwf_cleanup_threads();
#else //_SPREC
    fftw_cleanup_threads();
#endif //_SPREC
  }
#endif //_HAVE_OMP

#ifdef _SPREC
  fftwf_mpi_cleanup();
#else //_SPREC
  fftw_mpi_cleanup();
#endif //_SPREC

#else //_HAVE_MPI

#ifdef _HAVE_OMP
#ifdef _SPREC
  fftwf_cleanup_threads();
#else //_SPREC
  fftw_cleanup_threads();
#endif //_SPREC
#endif //_HAVE_OMP

#endif //_HAVE_MPI
}

static void create_density_velpot_and_npot_fourier(ParamCoLoRe *par)
{
  //////
  // Generates a random realization of the delta_k
  // from the linear P_k.

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,IThread0)
#endif //_HAVE_OMP
  {
    int ii;
    double dk=2*M_PI/par->l_box;
    double idk3=1./(dk*dk*dk);
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    unsigned int seed_thr=par->seed_rng+IThread0+ithr;
    gsl_rng *rng_thr=init_rng(seed_thr);
    double factor=par->fgrowth_0*par->hubble_0;
    double factor2=par->hubble_0*par->hubble_0*par->OmegaM;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<par->nz_here;ii++) {
      int jj,ii_true;
      double kz;
      ii_true=par->iz0_here+ii;
      if(2*ii_true<=par->n_grid)
	kz=ii_true*dk;
      else
	kz=-(par->n_grid-ii_true)*dk;
      for(jj=0;jj<par->n_grid;jj++) {
	int kk;
	double ky;
	if(2*jj<=par->n_grid)
	  ky=jj*dk;
	else
	  ky=-(par->n_grid-jj)*dk;
	for(kk=0;kk<=par->n_grid/2;kk++) {
	  double kx;
	  double k_mod2;
	  lint index=kk+(par->n_grid/2+1)*((lint)(jj+par->n_grid*ii)); //Grid index for +k
	  double delta_mod,delta_phase;
	  if(2*kk<=par->n_grid)
	    kx=kk*dk;
	  else
	    kx=-(par->n_grid-kk)*dk; //This should never happen

	  k_mod2=kx*kx+ky*ky+kz*kz;

	  if(k_mod2<=0) {
	    par->grid_dens_f[index]=0;
	    par->grid_vpot_f[index]=0;
      par->grid_npot_f[index]=0;
	  }
	  else {
	    double lgk=0.5*log10(k_mod2);
	    double sigma2=pk_linear0(par,lgk)*idk3;
	    if(par->do_smoothing)
	      sigma2*=exp(-par->r2_smooth*k_mod2);
	    rng_delta_gauss(&delta_mod,&delta_phase,rng_thr,sigma2);
	    par->grid_dens_f[index]=delta_mod*cexp(I*delta_phase);
	    par->grid_vpot_f[index]=par->grid_dens_f[index]*factor/k_mod2;
      //Newtonian potential Phi(k,z) -1/k^2 3/2 H0^2 OmegaM delta(k,z)
      par->grid_npot_f[index]=-1.5*par->grid_dens_f[index]*factor2/k_mod2;
	  }
	}
      }
    }
    end_rng(rng_thr);
  }
}

static void create_density_and_velpot_fourier(ParamCoLoRe *par)
{
  //////
  // Generates a random realization of the delta_k
  // from the linear P_k.

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,IThread0)
#endif //_HAVE_OMP
  {
    int ii;
    double dk=2*M_PI/par->l_box;
    double idk3=1./(dk*dk*dk);
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    unsigned int seed_thr=par->seed_rng+IThread0+ithr;
    gsl_rng *rng_thr=init_rng(seed_thr);
    double factor=par->fgrowth_0*par->hubble_0;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<par->nz_here;ii++) {
      int jj,ii_true;
      double kz;
      ii_true=par->iz0_here+ii;
      if(2*ii_true<=par->n_grid)
	kz=ii_true*dk;
      else
	kz=-(par->n_grid-ii_true)*dk;
      for(jj=0;jj<par->n_grid;jj++) {
	int kk;
	double ky;
	if(2*jj<=par->n_grid)
	  ky=jj*dk;
	else
	  ky=-(par->n_grid-jj)*dk;
	for(kk=0;kk<=par->n_grid/2;kk++) {
	  double kx;
	  double k_mod2;
	  lint index=kk+(par->n_grid/2+1)*((lint)(jj+par->n_grid*ii)); //Grid index for +k
	  double delta_mod,delta_phase;
	  if(2*kk<=par->n_grid)
	    kx=kk*dk;
	  else
	    kx=-(par->n_grid-kk)*dk; //This should never happen

	  k_mod2=kx*kx+ky*ky+kz*kz;

	  if(k_mod2<=0) {
	    par->grid_dens_f[index]=0;
	    par->grid_vpot_f[index]=0;
	  }
	  else {
	    double lgk=0.5*log10(k_mod2);
	    double sigma2=pk_linear0(par,lgk)*idk3;
	    if(par->do_smoothing)
	      sigma2*=exp(-par->r2_smooth*k_mod2);
	    rng_delta_gauss(&delta_mod,&delta_phase,rng_thr,sigma2);
	    par->grid_dens_f[index]=delta_mod*cexp(I*delta_phase);
	    par->grid_vpot_f[index]=par->grid_dens_f[index]*factor/k_mod2;
  	  }
	}
      }
    }
    end_rng(rng_thr);
  }
}

static void radial_velocity_from_potential(ParamCoLoRe *par)
{
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par)
#endif //_HAVE_OMP
  {
    double dx=par->l_box/par->n_grid;
    double idx=1./dx;
    int iz;
    int ngx=2*(par->n_grid/2+1);

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint iz_hi=iz+1;
      lint iz_lo=iz-1;
      lint iz_0=iz;
      double z=dx*(iz+par->iz0_here+0.5)-par->pos_obs[2];
      if(iz==0) iz_lo=par->n_grid-1;
      if(iz==par->n_grid-1) iz_hi=0;
      iz_hi*=ngx*par->n_grid;
      iz_lo*=ngx*par->n_grid;
      iz_0*=ngx*par->n_grid;
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint iy_hi=iy+1;
	lint iy_lo=iy-1;
	lint iy_0=iy;
	double y=dx*(iy+0.5)-par->pos_obs[1];
	if(iy==0) iy_lo=par->n_grid-1;
	if(iy==par->n_grid-1) iy_hi=0;
	iy_hi*=ngx;
	iy_lo*=ngx;
	iy_0*=ngx;
	for(ix=0;ix<par->n_grid;ix++) {
	  double vel[3];
	  double ur[3];
	  lint ix_hi=ix+1;
	  lint ix_lo=ix-1;
	  lint ix_0=ix;
	  double x=dx*(ix+0.5)-par->pos_obs[0];
	  double irr=1./sqrt(x*x+y*y+z*z);
	  if(ix==0) ix_lo=par->n_grid-1;
	  if(ix==par->n_grid-1) ix_hi=0;

	  ur[0]=x*irr;
	  ur[1]=y*irr;
	  ur[2]=z*irr;

	  vel[0]=0.5*idx*(par->grid_vpot[ix_hi+iy_0+iz_0]-par->grid_vpot[ix_lo+iy_0+iz_0]);
	  vel[1]=0.5*idx*(par->grid_vpot[ix_0+iy_hi+iz_0]-par->grid_vpot[ix_0+iy_lo+iz_0]);
	  if(iz==0)
	    vel[2]=0.5*idx*(par->grid_vpot[ix_0+iy_0+iz_hi]-par->slice_left[ix_0+iy_0]);
	  else if(iz==par->nz_here-1)
	    vel[2]=0.5*idx*(par->slice_right[ix_0+iy_0]-par->grid_vpot[ix_0+iy_0+iz_lo]);
	  else
	    vel[2]=0.5*idx*(par->grid_vpot[ix_0+iy_0+iz_hi]-par->grid_vpot[ix_0+iy_0+iz_lo]);

	  par->grid_rvel[ix_0+iy_0+iz_0]=vel[0]*ur[0]+vel[1]*ur[1]+vel[2]*ur[2];
	}
      }
    } // end omp for
  } // end omp parallel
}

void create_d_and_vr_fields(ParamCoLoRe *par)
{
  //////
  // Creates a realization of the gaussian density
  // contrast field from the linear P_k

  lint n_grid_tot=2*(par->n_grid/2+1)*((lint)(par->n_grid*par->nz_here));
  print_info("*** Creating Gaussian density field \n");

  print_info("Creating Fourier-space density and velocity potential \n");
  if(NodeThis==0) timer(0);
  create_density_and_velpot_fourier(par);
  if(NodeThis==0) timer(2);

  print_info("Transforming density and velocity potential\n");
  if(NodeThis==0) timer(0);
  fftw_wrap(par->n_grid,par->grid_dens_f,par->grid_dens);
  fftw_wrap(par->n_grid,par->grid_vpot_f,par->grid_vpot);
  if(NodeThis==0) timer(2);

  print_info("Normalizing density and velocity potential \n");
  if(NodeThis==0) timer(0);
#ifdef _HAVE_OMP
#pragma omp parallel default(none)				\
  shared(n_grid_tot,par)
#endif //_HAVE_OMP
  {
    lint ii;
    double norm=pow(sqrt(2*M_PI)/par->l_box,3);

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<n_grid_tot;ii++) {
      par->grid_dens[ii]*=norm;
      par->grid_vpot[ii]*=norm;
    }
  }//end omp parallel
  if(NodeThis==0) timer(2);

  lint slice_size=2*(par->n_grid/2+1)*par->n_grid;
#ifdef _HAVE_MPI
  MPI_Status stat;

  //Pass rightmost slice to right node and receive left slice from left node
  MPI_Sendrecv(&(par->grid_vpot[(par->nz_here-1)*slice_size]),slice_size,FLOUBLE_MPI,NodeRight,1,
	       par->slice_left,slice_size,FLOUBLE_MPI,NodeLeft,1,MPI_COMM_WORLD,&stat);
  //Pass leftmost slice to left node and receive right slice from right node
  MPI_Sendrecv(par->grid_vpot,slice_size,FLOUBLE_MPI,NodeLeft,2,
	       par->slice_right,slice_size,FLOUBLE_MPI,NodeRight,2,MPI_COMM_WORLD,&stat);
#else //_HAVE_MPI
  par->slice_left=&(par->grid_vpot[(par->n_grid-1)*slice_size]);
  par->slice_right=par->grid_vpot;
#endif //_HAVE_MPI

  print_info("Calculating radial velocity \n");
  if(NodeThis==0) timer(0);
  radial_velocity_from_potential(par);
  if(NodeThis==0) timer(2);

  compute_sigma_dens(par);

  print_info("\n");

  //Output density field if necessary
  if(par->output_density)
    write_grid(par);
}
void create_d_phi_and_vr_fields(ParamCoLoRe *par)
{
  //////
  // Creates a realization of the gaussian density
  // contrast field from the linear P_k

  lint n_grid_tot=2*(par->n_grid/2+1)*((lint)(par->n_grid*par->nz_here));
  print_info("*** Creating Gaussian density field \n");

  print_info("Creating Fourier-space density and velocity potential \n");
  if(NodeThis==0) timer(0);
  create_density_velpot_and_npot_fourier(par);
  if(NodeThis==0) timer(2);

  print_info("Transforming density and velocity potential\n");
  if(NodeThis==0) timer(0);
  fftw_wrap(par->n_grid,par->grid_dens_f,par->grid_dens);
  fftw_wrap(par->n_grid,par->grid_vpot_f,par->grid_vpot);
  fftw_wrap(par->n_grid,par->grid_npot_f,par->grid_npot);
  if(NodeThis==0) timer(2);

  print_info("Normalizing density and velocity potential \n");
  if(NodeThis==0) timer(0);
#ifdef _HAVE_OMP
#pragma omp parallel default(none)				\
  shared(n_grid_tot,par)
#endif //_HAVE_OMP
  {
    lint ii;
    double norm=pow(sqrt(2*M_PI)/par->l_box,3);

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<n_grid_tot;ii++) {
      par->grid_dens[ii]*=norm;
      par->grid_vpot[ii]*=norm;
      par->grid_npot[ii]*=norm;
    }
  }//end omp parallel
  if(NodeThis==0) timer(2);

  lint slice_size=2*(par->n_grid/2+1)*par->n_grid;
#ifdef _HAVE_MPI
  MPI_Status stat;

  //Pass rightmost slice to right node and receive left slice from left node
  MPI_Sendrecv(&(par->grid_vpot[(par->nz_here-1)*slice_size]),slice_size,FLOUBLE_MPI,NodeRight,1,
	       par->slice_left,slice_size,FLOUBLE_MPI,NodeLeft,1,MPI_COMM_WORLD,&stat);
  //Pass leftmost slice to left node and receive right slice from right node
  MPI_Sendrecv(par->grid_vpot,slice_size,FLOUBLE_MPI,NodeLeft,2,
	       par->slice_right,slice_size,FLOUBLE_MPI,NodeRight,2,MPI_COMM_WORLD,&stat);
#else //_HAVE_MPI
  par->slice_left=&(par->grid_vpot[(par->n_grid-1)*slice_size]);
  par->slice_right=par->grid_vpot;
#endif //_HAVE_MPI

  print_info("Calculating radial velocity \n");
  if(NodeThis==0) timer(0);
  radial_velocity_from_potential(par);
  if(NodeThis==0) timer(2);

  compute_sigma_dens(par);

  print_info("\n");

  //Output density field if necessary
  if(par->output_density)
    write_grid(par);
  if(par->output_potential)
    write_pot(par);
}
