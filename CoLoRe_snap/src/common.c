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

typedef struct {
  int i;
  flouble d;
} IsortStruct;

static int compareIsort(const void *a,const void *b)
{
  IsortStruct *ia=(IsortStruct *)a;
  IsortStruct *ib=(IsortStruct *)b;
  if(ia->d<ib->d) return -1;
  else return 1;
}

int *ind_sort(int n,flouble *arr)
{
  int i;
  IsortStruct *st=my_malloc(n*sizeof(IsortStruct));
  for(i=0;i<n;i++) {
    st[i].i=i;
    st[i].d=arr[i];
  }

  qsort(st,n,sizeof(IsortStruct),compareIsort);

  int *iarr=my_malloc(n*sizeof(int));
  for(i=0;i<n;i++)
    iarr[i]=st[i].i;

  free(st);
  return iarr;
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

void rng_gauss(gsl_rng *rng,double *r1,double *r2)
{
  double phase=2*M_PI*rng_01(rng);
  double u=sqrt(-2*log(1-rng_01(rng)));
  *r1=u*cos(phase);
  *r2=u*sin(phase);
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
  MPIThreadsOK=provided>=MPI_THREAD_FUNNELED;
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
      printf("MPI task %d has %d OMP threads\n",ii,nthreads_all[ii]);
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
  printf("MPI task %d, OMP thread count starts at %d\n",NodeThis,IThread0);
  print_info(" MPIThreadsOK = %d\n",MPIThreadsOK);
#endif //_DEBUG
}

CatalogCartesian *catalog_cartesian_alloc(int nsrcs)
{
  CatalogCartesian *cat=my_malloc(sizeof(CatalogCartesian));

  if(nsrcs>0) {
    cat->nsrc=nsrcs;
    cat->pos=my_malloc(NPOS_CC*nsrcs*sizeof(float));
  }
  else {
    cat->nsrc=0;
    cat->pos=NULL;
  }

  return cat;
}

void catalog_cartesian_free(CatalogCartesian *cat)
{
  if(cat->nsrc>0) {
    free(cat->pos);
  }
  free(cat);
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
    fprintf(stderr,"MPI task %d, Fatal: %s",NodeThis,msg);
    exit(level);
  }
  else {
    fprintf(stderr,"MPI task %d, Warning: %s",NodeThis,msg);
  }
}

size_t my_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  if(fwrite(ptr,size,nmemb,stream)!=nmemb)
    report_error(1,"Error fwriting\n");

  return nmemb;
}

unsigned long long get_max_memory(ParamCoLoRe *par,int just_test)
{
  unsigned long long total_GB=0;
  unsigned long long total_GB_gau=0;
  unsigned long long total_GB_lpt=0;
  int fac_gau=2;

  total_GB_gau=(fac_gau*(par->nz_here+1)*((long)(par->n_grid*(par->n_grid/2+1))))*sizeof(dftw_complex);

  if(par->dens_type==DENS_TYPE_1LPT) {
    total_GB_lpt=(unsigned long long)(3*(1+par->lpt_buffer_fraction)*par->nz_here*
				      ((long)((par->n_grid/2+1)*par->n_grid))*sizeof(dftw_complex));
  }
  else if(par->dens_type==DENS_TYPE_2LPT) {
    total_GB_lpt=0;
    total_GB_lpt=(unsigned long long)(8*(1+par->lpt_buffer_fraction)*par->nz_here*
				      ((long)((par->n_grid/2+1)*par->n_grid))*sizeof(dftw_complex));
  }

  unsigned long long total_GB_srcs=0;
  if(par->do_srcs) {
    int ipop;
    long nsrc=0;
    double vol=par->l_box*par->l_box*par->l_box;
    for(ipop=0;ipop<par->n_srcs;ipop++) {
      long ngal=(long)(par->ndens[ipop]*vol);
      nsrc+=ngal;
    }
    long size_source=NPOS_CC*sizeof(flouble)+sizeof(int);
    total_GB_srcs=size_source*nsrc;
  }

  total_GB=total_GB_gau+total_GB_lpt+total_GB_srcs;

#ifdef _DEBUG
  int jj;
  for(jj=0;jj<NNodes;jj++) {
    if(jj==NodeThis) {
      printf("Node %d will allocate %.3lf GB [",NodeThis,(double)(total_GB/pow(1024.,3)));
      printf("%.3lf GB (Gaussian)",(double)(total_GB_gau/pow(1024.,3)));
      if((par->dens_type==DENS_TYPE_1LPT) || (par->dens_type==DENS_TYPE_2LPT))
	printf(", %.3lf GB (%dLPT)",(double)(total_GB_lpt/pow(1024.,3)),par->dens_type);
      if(par->do_srcs)
	printf(", %.3lf GB (srcs)",(double)(total_GB_srcs/pow(1024.,3)));
      printf("]\n");
    }
#ifdef _HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif //_HAVE_MPI
  }
#endif //_DEBUG

  if(just_test==0) {
    void *ptest=my_malloc(total_GB);
    free(ptest);
  }

  return total_GB;
}
