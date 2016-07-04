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

#ifdef _HAVE_OMP
static double relbeg,relend,absbeg,absend;
#else //_HAVE_OMP
static time_t relbeg,relend,absbeg,absend;
#endif //_HAVE_OMP
int NodeThis=0;
int NodeLeft=0;
int NodeRight=0;
int NNodes=1;
int IThread0=0;
int MPIThreadsOK=0;

void *my_malloc(size_t size)
{
  void *ptrout=malloc(size);
  if(ptrout==NULL) {
    fprintf(stderr,"out of memory\n");
    exit(1);
  }
  
  return ptrout;
}

void *my_calloc(size_t nmemb,size_t size)
{
  void *ptrout=calloc(nmemb,size);
  if(ptrout==NULL) {
    fprintf(stderr,"out of memory\n");
    exit(1);
  }

  return ptrout;
}

void error_open_file(char *fname)
{
  fprintf(stderr,"CoLoRe: Couldn't open file %s \n",fname);
  exit(1);
}

void error_read_line(char *fname,int nlin)
{
  fprintf(stderr,"CoLoRe: Error reading file %s, line %d \n",fname,nlin);
  exit(1);
}

int linecount(FILE *f)
{
  //////
  // Counts #lines from file
  int i0=0;
  char ch[1000];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

void timer(int i)
{
  /////
  // Timing routine
  // timer(0) -> initialize relative clock
  // timer(1) -> read relative clock
  // timer(2) -> read relative clock and initialize it afterwards
  // timer(4) -> initialize absolute clock
  // timer(5) -> read absolute clock
#ifdef _HAVE_OMP
  if(i==0)
    relbeg=omp_get_wtime();
  else if(i==1) {
    relend=omp_get_wtime();
    printf(">    Relative time ellapsed %.1lf ms\n",1000*(relend-relbeg));
  }
  else if(i==2) {
    relend=omp_get_wtime();
    printf(">    Relative time ellapsed %.1lf ms\n",1000*(relend-relbeg));
    relbeg=omp_get_wtime();
  }
  else if(i==4)
    absbeg=omp_get_wtime();
  else if(i==5) {
    absend=omp_get_wtime();
    printf(">    Total time ellapsed %.1lf ms\n",1000*(absend-absbeg));
  }
#else //_NO_OMP
  int diff;
  
  if(i==0)
    relbeg=time(NULL);
  else if(i==1) {
    relend=time(NULL);
    diff=(int)(difftime(relend,relbeg));
    printf(">    Relative time ellapsed %02d:%02d:%02d \n",
	   diff/3600,(diff/60)%60,diff%60);
  }
  else if(i==2) {
    relend=time(NULL);
    diff=(int)(difftime(relend,relbeg));
    printf(">    Relative time ellapsed %02d:%02d:%02d \n",
	   diff/3600,(diff/60)%60,diff%60);
    relbeg=time(NULL);
  }
  else if(i==4)
    absbeg=time(NULL);
  else if(i==5) {
    absend=time(NULL);
    diff=(int)(difftime(absend,absbeg));
    printf(">    Total time ellapsed %02d:%02d:%02d \n",
	   diff/3600,(diff/60)%60,diff%60);
  }
#endif //_NO_OMP
}

gsl_rng *init_rng(unsigned int seed)
{
  //  gsl_rng *rng=gsl_rng_alloc(gsl_rng_ranlux);
  gsl_rng *rng=gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng,seed);

  return rng;
}

double rng_01(gsl_rng *rng)
{
  double result=gsl_rng_uniform(rng);
  return result;
}

int rng_poisson(double lambda,gsl_rng *rng)
{
  unsigned int pois=gsl_ran_poisson(rng,lambda);
  return (int)pois;
}

void rng_delta_gauss(double *module,double *phase,
		     gsl_rng *rng,double sigma2)
{
  //////
  // Returns module and phase of two random 
  // gaussian numbers. I.e.: 
  double u;
  *phase=2*M_PI*rng_01(rng);
  u=rng_01(rng);
  *module=sqrt(-sigma2*log(1-u));
}

void end_rng(gsl_rng *rng)
{
  gsl_rng_free(rng);
}

void mpi_init(int* p_argc,char*** p_argv)
{
#ifdef _HAVE_MPI
  int ii,nthreads_this;
  int *nthreads_all;
#ifdef _HAVE_OMP
  int provided;
  MPI_Init_thread(p_argc,p_argv,MPI_THREAD_FUNNELED,&provided);
  MPIThreadsOK = provided >= MPI_THREAD_FUNNELED;
#else //_HAVE_OMP
  MPI_Init(p_argc,p_argv);
  MPIThreadsOK=0;
#endif //_HAVE_OMP

  MPI_Comm_size(MPI_COMM_WORLD,&NNodes);
  MPI_Comm_rank(MPI_COMM_WORLD,&NodeThis);
  if(NodeThis==0)
    NodeLeft=NNodes-1;
  else
    NodeLeft=NodeThis-1;
  if(NodeThis==NNodes-1)
    NodeRight=0;
  else
    NodeRight=NodeThis+1;

  nthreads_all=my_malloc(NNodes*sizeof(int));
#ifdef _HAVE_OMP
  nthreads_this=omp_get_max_threads();
#else //_HAVE_OMP
  nthreads_this=1;
#endif //_HAVE_OMP
  MPI_Allgather(&nthreads_this,1,MPI_INT,nthreads_all,1,MPI_INT,MPI_COMM_WORLD);
#ifdef _DEBUG
  if(NodeThis==0) {
    for(ii=0;ii<NNodes;ii++)
      printf("Node %d has %d threads\n",ii,nthreads_all[ii]);
  }
#endif //_DEBUG
  IThread0=0;
  for(ii=0;ii<NodeThis;ii++)
    IThread0+=nthreads_all[ii];
  free(nthreads_all);

#else //_HAVE_MPI
  NodeThis=0;
  NNodes=1;
  IThread0=0;
#ifdef _HAVE_OMP
  MPIThreadsOK=1;
#else //_HAVE_OMP
  MPIThreadsOK=0;
#endif //_HAVE_OMP
#endif //_HAVE_MPI
#ifdef _DEBUG
  printf("Node %d, thread count starts at %d\n",NodeThis,IThread0);
  print_info(" Threads = %d\n",MPIThreadsOK);
#endif //_DEBUG
}

void print_info(char *fmt,...)
{
  if(NodeThis==0) {
    va_list args;
    char msg[256];
    
    va_start(args,fmt);
    vsprintf(msg,fmt,args);
    va_end(args);
    
    printf("%s",msg);
  }
}

void report_error(int level,char *fmt,...)
{
  va_list args;
  char msg[256];

  va_start(args,fmt);
  vsprintf(msg,fmt,args);
  va_end(args);
  
  if(level) {
    fprintf(stderr,"Node %d, Fatal: %s",NodeThis,msg);
    exit(level);
  }
  else {
    fprintf(stderr,"Node %d, Warning: %s",NodeThis,msg);
  }
}

size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  if(fwrite(ptr,size,nmemb,stream)!=nmemb)
    report_error(1,"Error fwriting\n");

  return nmemb;
}

#ifndef _OLD_IM
void alloc_onion_info(ParamCoLoRe *par,OnionInfo *oi,
		      int nside,int nz,
		      flouble *z0_arr,flouble *zf_arr,
		      int add_rsd,int add_next)
{
  int izz;
  long npix_ang=nside2npix(nside);
  double dx=par->l_box/par->n_grid;
  double dz_additional,z0_here,zf_here;
  double dr_perp_max=0,dr_par_max=0;

  oi->r0_arr=my_malloc(nz*sizeof(flouble));
  oi->rf_arr=my_malloc(nz*sizeof(flouble));

  for(izz=0;izz<nz;izz++) {
    double r0=r_of_z(par,z0_arr[izz]);
    double rf=r_of_z(par,zf_arr[izz]);
    double dr_par=(rf-r0);
    double dr_perp=(rf+r0)/2*sqrt(4*M_PI/npix_ang);
    oi->r0_arr[izz]=r0;
    oi->rf_arr[izz]=rf;
    if(dr_par>dr_par_max)
      dr_par_max=dr_par;
    if(dr_perp>dr_perp_max)
      dr_perp_max=dr_perp;
#ifdef _DEBUG
    print_info("%d z [%lf,%lf], r=%lf, dr_par=%lf, dr_perp=%lf\n",
	       izz,z0_arr[izz],zf_arr[izz],(r0+rf)/2,dr_par,dr_perp);
#endif //_DEBUG
  }
  dz_additional=dr_perp_max*dr_perp_max;
  if(add_rsd)
    dz_additional+=DR_RSD_ADDITIONAL*DR_RSD_ADDITIONAL;
  dz_additional=sqrt(dz_additional);
  if(add_next)
    dz_additional+=dx;

#ifdef _DEBUG
  print_info("Will use a buffer of %lf Mpc/h\n",dz_additional);
#endif //_DEBUG
  z0_here=dx*par->iz0_here-par->pos_obs[2]-dz_additional;
  zf_here=dx*(par->nz_here+par->iz0_here)-par->pos_obs[2]+dz_additional;
  oi->list_ipix=my_malloc(nz*sizeof(long *));
  oi->num_pix=my_malloc(nz*sizeof(int));
  oi->maps=my_malloc(nz*sizeof(flouble *));
  oi->nz=nz;

  //TODO: omp this
#pragma omp parallel default(none)			\
  shared(par,z0_arr,zf_arr,oi,npix_ang,nz,nside)	\
  shared(z0_here,zf_here)
  {
    int iz;
    
#pragma omp for
    for(iz=0;iz<nz;iz++) {
      long ipx;
      double r0=oi->r0_arr[iz];
      double rf=oi->rf_arr[iz];
      oi->num_pix[iz]=0;
      for(ipx=0;ipx<npix_ang;ipx++) {
	double pos[3];
	double zpix0,zpixf;
	//TODO: is this better than pix2ang?
	pix2vec_ring(nside,ipx,pos);
	zpix0=r0*pos[2];
	zpixf=rf*pos[2];
	if(((zpix0>=z0_here) && (zpix0<=zf_here)) ||
	   ((zpixf>=z0_here) && (zpixf<=zf_here))) {
	  oi->num_pix[iz]++;
	}
      }
    } //end omp for
  } //end omp parallel

  for(izz=0;izz<nz;izz++) {
    if(oi->num_pix[izz]>0) {
      oi->list_ipix[izz]=my_malloc(oi->num_pix[izz]*sizeof(long));
      oi->maps[izz]=my_calloc(oi->num_pix[izz],sizeof(flouble));
    }
  }

#pragma omp parallel default(none)			\
  shared(par,z0_arr,zf_arr,oi,npix_ang,nz,nside)	\
  shared(z0_here,zf_here)
  {
    int iz;
    
#pragma omp for
    for(iz=0;iz<nz;iz++) {
      long ipx;
      double r0=oi->r0_arr[iz];
      double rf=oi->rf_arr[iz];
      oi->num_pix[iz]=0;
      for(ipx=0;ipx<npix_ang;ipx++) {
	double pos[3];
	double zpix0,zpixf;
	pix2vec_ring(nside,ipx,pos);
	zpix0=r0*pos[2];
	zpixf=rf*pos[2];
	if(((zpix0>=z0_here) && (zpix0<=zf_here)) ||
	   ((zpixf>=z0_here) && (zpixf<=zf_here))) {
	  oi->list_ipix[iz][oi->num_pix[iz]]=ipx;
	  oi->num_pix[iz]++;
	}
      }
    } //end omp for
  } //end omp parallel
}

void free_onion_info(OnionInfo *oi)
{
  int iz;
  for(iz=0;iz<oi->nz;iz++) {
    if(oi->num_pix[iz]>0) {
      free(oi->maps[iz]);
      free(oi->list_ipix[iz]);
    }
  }
  free(oi->r0_arr);
  free(oi->rf_arr);
  free(oi->maps);
  free(oi->list_ipix);
  free(oi->num_pix);
}
#endif //_OLD_IM
