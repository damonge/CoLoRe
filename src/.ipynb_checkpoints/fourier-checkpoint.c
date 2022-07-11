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
  //  double s2_save=par->sigma2_gauss;
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
    long ng_tot=par->n_grid*((long)(par->n_grid*par->n_grid));

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long iz0=iz*((long)(2*(par->n_grid/2+1)*par->n_grid));
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long iy0=iy*2*(par->n_grid/2+1);
	for(ix=0;ix<par->n_grid;ix++) {
	  long index=ix+iy0+iz0;
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
  //  if(par->need_onions)
  //    par->sigma2_gauss=s2_save;
}

void fftw_wrap_c2r(int ng,dftw_complex *pin,flouble *pout)
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

void fftw_wrap_r2c(int ng,flouble *pin,dftw_complex *pout)
{
#ifdef _SPREC
  fftwf_plan plan_ft;
#ifdef _HAVE_MPI
  plan_ft=fftwf_mpi_plan_dft_r2c_3d(ng,ng,ng,pin,pout,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else //_HAVE_MPI
  plan_ft=fftwf_plan_dft_r2c_3d(ng,ng,ng,pin,pout,FFTW_ESTIMATE);
#endif //_HAVE_MPI
  fftwf_execute(plan_ft);
  fftwf_destroy_plan(plan_ft);
#else //_SPREC
  fftw_plan plan_ft;
#ifdef _HAVE_MPI
  plan_ft=fftw_mpi_plan_dft_r2c_3d(ng,ng,ng,pin,pout,MPI_COMM_WORLD,FFTW_ESTIMATE);
#else //_HAVE_MPI
  plan_ft=fftw_plan_dft_r2c_3d(ng,ng,ng,pin,pout,FFTW_ESTIMATE);
#endif //_HAVE_MPI
  fftw_execute(plan_ft);
  fftw_destroy_plan(plan_ft);
#endif //_SPREC
}

void init_fftw(ParamCoLoRe *par)
{
  //Set FFTW domain decomposition
  par->nz_all =my_calloc(NNodes,sizeof(int));
  par->iz0_all=my_calloc(NNodes,sizeof(int));

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
  fftwf_mpi_local_size_3d(par->n_grid,par->n_grid,par->n_grid/2+1,MPI_COMM_WORLD,&nz,&iz0);
#else //_SPREC
  fftw_mpi_local_size_3d(par->n_grid,par->n_grid,par->n_grid/2+1,MPI_COMM_WORLD,&nz,&iz0);
#endif //_SPREC
  par->nz_here=nz;
  par->iz0_here=iz0;

  MPI_Allreduce(&(par->nz_here),&(par->nz_max),1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allgather(&(par->nz_here),1,MPI_INT,par->nz_all,1,MPI_INT,MPI_COMM_WORLD);
  MPI_Allgather(&(par->iz0_here),1,MPI_INT,par->iz0_all,1,MPI_INT,MPI_COMM_WORLD);

#else //_HAVE_MPI

#ifdef _HAVE_OMP
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
#ifdef _SPREC
  fftwf_plan_with_nthreads(omp_get_max_threads());
#else //_SPREC
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif //_SPREC
#endif //_HAVE_OMP

  par->nz_here=par->n_grid;
  par->iz0_here=0;
  par->nz_max=par->nz_here;
  par->nz_all[0]=par->nz_here;
  par->iz0_all[0]=par->iz0_here;
#endif //_HAVE_MPI
}

void allocate_fftw(ParamCoLoRe *par)
{
  //Allocate all memory for grids
  ptrdiff_t dsize=par->nz_max*((long)(par->n_grid*(par->n_grid/2+1)));

#ifdef _SPREC
  par->grid_dens_f=fftwf_alloc_complex(dsize);
#else //_SPREC
  par->grid_dens_f=fftw_alloc_complex(dsize);
#endif //_SPREC
  if(par->grid_dens_f==NULL)
    report_error(1,"Ran out of memory\n");
  par->grid_dens=(flouble *)(par->grid_dens_f);

#ifdef _SPREC
  par->grid_npot_f=fftwf_alloc_complex(dsize+2*par->n_grid*(par->n_grid/2+1));
#else //_SPREC
  par->grid_npot_f= fftw_alloc_complex(dsize+2*par->n_grid*(par->n_grid/2+1));
#endif //_SPREC
  if(par->grid_npot_f==NULL)
    report_error(1,"Ran out of memory\n");
  par->grid_npot=(flouble *)(par->grid_npot_f);

#ifdef _HAVE_MPI
  par->slice_left =&(par->grid_npot[2*dsize]);
  par->slice_right=&(par->grid_npot[2*(dsize+par->n_grid*(par->n_grid/2+1))]);
#endif //_HAVE_MPI
}

void end_fftw(ParamCoLoRe *par)
{
#ifdef _SPREC
  if(par->grid_dens_f!=NULL)
    fftwf_free(par->grid_dens_f);
  if(par->grid_npot_f!=NULL)
    fftwf_free(par->grid_npot_f);
#else //_SPREC
  if(par->grid_dens_f!=NULL)
    fftw_free(par->grid_dens_f);
  if(par->grid_npot_f!=NULL)
    fftw_free(par->grid_npot_f);
#endif //_SPREC

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

static void create_grids_fourier(ParamCoLoRe *par)
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
	  long index=kk+(par->n_grid/2+1)*((long)(jj+par->n_grid*ii)); //Grid index for +k
	  double delta_mod,delta_phase;
	  if(2*kk<=par->n_grid)
	    kx=kk*dk;
	  else
	    kx=-(par->n_grid-kk)*dk; //This should never happen
	  
	  k_mod2=kx*kx+ky*ky+kz*kz;
	  
	  if(k_mod2<=0) {
	    par->grid_dens_f[index]=0;
	    par->grid_npot_f[index]=0;
	  }
	  else {
	    double lgk=0.5*log10(k_mod2);
	    double sigma2=pk_linear0(par,lgk)*idk3;
	    rng_delta_gauss(&delta_mod,&delta_phase,rng_thr,sigma2);
	    par->grid_dens_f[index]=delta_mod*cexp(I*delta_phase);
	    par->grid_npot_f[index]=-par->prefac_lensing*par->grid_dens_f[index]/k_mod2;
	    if(par->do_smoothing) {
	      double sm=exp(-0.5*par->r2_smooth*k_mod2);
	      par->grid_dens_f[index]*=sm;
	      if(par->smooth_potential)
		par->grid_npot_f[index]*=sm;
	    }
	  }
	}
      }
    }
    end_rng(rng_thr);
  }
}

void create_cartesian_fields(ParamCoLoRe *par)
{
  //////
  // Creates a realization of the gaussian density
  // contrast field from the linear P_k
  
  long n_grid_tot=2*(par->n_grid/2+1)*((long)(par->n_grid*par->nz_here));
  print_info("*** Creating Gaussian density field \n");

  print_info("Creating Fourier-space density and Newtonian potential \n");
  if(NodeThis==0) timer(0);
  create_grids_fourier(par);
  if(NodeThis==0) timer(2);

  print_info("Transforming density and Newtonian potential\n");
  if(NodeThis==0) timer(0);
  fftw_wrap_c2r(par->n_grid,par->grid_dens_f,par->grid_dens);
  fftw_wrap_c2r(par->n_grid,par->grid_npot_f,par->grid_npot);
  if(NodeThis==0) timer(2);

  print_info("Normalizing density and Newtonian potential \n");
  if(NodeThis==0) timer(0);
#ifdef _HAVE_OMP
#pragma omp parallel default(none)				\
  shared(n_grid_tot,par)
#endif //_HAVE_OMP
  {
    long ii;
    double norm=pow(sqrt(2*M_PI)/par->l_box,3);
    
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<n_grid_tot;ii++) {
      par->grid_dens[ii]*=norm;
      par->grid_npot[ii]*=norm;
    }
  }//end omp parallel
  if(NodeThis==0) timer(2);

  long slice_size=2*(par->n_grid/2+1)*par->n_grid;
#ifdef _HAVE_MPI
  MPI_Status stat;

  //Pass rightmost slice to right node and receive left slice from left node
  MPI_Sendrecv(&(par->grid_npot[(par->nz_here-1)*slice_size]),slice_size,FLOUBLE_MPI,NodeRight,1,
	       par->slice_left,slice_size,FLOUBLE_MPI,NodeLeft,1,MPI_COMM_WORLD,&stat);
  //Pass leftmost slice to left node and receive right slice from right node
  MPI_Sendrecv(par->grid_npot,slice_size,FLOUBLE_MPI,NodeLeft,2,
	       par->slice_right,slice_size,FLOUBLE_MPI,NodeRight,2,MPI_COMM_WORLD,&stat);
#else //_HAVE_MPI
  par->slice_left=&(par->grid_npot[(par->n_grid-1)*slice_size]);
  par->slice_right=par->grid_npot;
#endif //_HAVE_MPI

  compute_sigma_dens(par);

  print_info("\n");

  //Output density field if necessary
  if(par->output_density)
    write_density_grid(par,"gaussian");
}
