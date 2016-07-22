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
  printf("Node %d has %d threads\n",NodeThis,nthreads_all[NodeThis]);
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

OnionInfo *alloc_onion_info(ParamCoLoRe *par,int nside_base,flouble dr)
{
  int irr;
  double dx=par->l_box/par->n_grid;
  double z0_node=dx*par->iz0_here-par->pos_obs[2];
  double zf_node=dx*(par->iz0_here+par->nz_here)-par->pos_obs[2];
  OnionInfo *oi=my_malloc(sizeof(OnionInfo));

  oi->nr=(int)(par->r_max/dr)+1;
  oi->r0_arr=my_malloc(oi->nr*sizeof(flouble));
  oi->rf_arr=my_malloc(oi->nr*sizeof(flouble));
  oi->nside_arr=my_malloc(oi->nr*sizeof(int));
  oi->num_pix=my_malloc(oi->nr*sizeof(int));
  oi->list_ipix=my_malloc(oi->nr*sizeof(long *));
  oi->maps=my_malloc(oi->nr*sizeof(flouble *));

  for(irr=0;irr<oi->nr;irr++) {
    int nside_here=nside_base;
    flouble rm=(irr+0.5)*dr;
    flouble dr_trans=rm*sqrt(4*M_PI/(12*nside_here*nside_here));
    oi->r0_arr[irr]=irr*dr;
    oi->rf_arr[irr]=(irr+1)*dr;

    while(dr_trans>1.333*dr) {
      nside_here*=2;
      dr_trans=rm*sqrt(4*M_PI/(12*nside_here*nside_here));
    }
    oi->nside_arr[irr]=nside_here;
  }

#pragma omp parallel default(none)			\
  shared(oi,z0_node,zf_node,dx)
  {
    int ir;

#pragma omp for schedule(dynamic,DYNAMIC_SIZE)
    for(ir=0;ir<oi->nr;ir++) {
      double r0=oi->r0_arr[ir];
      double rf=oi->rf_arr[ir];
      double rm=0.5*(r0+rf);
      long ipx,npix=12*oi->nside_arr[ir]*oi->nside_arr[ir];
      flouble dr_trans=rm*sqrt(4*M_PI/npix);
      double dz_additional=dr_trans+dx;
      double z0_here=z0_node-dz_additional;
      double zf_here=zf_node+dz_additional;
      oi->num_pix[ir]=0;

      for(ipx=0;ipx<npix;ipx++) {
	double zpix0,zpixf,pos[3];
	//TODO: is this better than pix2ang?
	pix2vec_ring(oi->nside_arr[ir],ipx,pos);
	zpix0=r0*pos[2];
	zpixf=rf*pos[2];
	if(((zpix0>=z0_here) && (zpix0<=zf_here)) ||
	   ((zpixf>=z0_here) && (zpixf<=zf_here))) {
	  oi->num_pix[ir]++;
	}
      }
    } //end omp for
  } //end omp parallel

  for(irr=0;irr<oi->nr;irr++) {
    if(oi->num_pix[irr]>0) {
      oi->list_ipix[irr]=my_malloc(oi->num_pix[irr]*sizeof(long));
      oi->maps[irr]=my_calloc(oi->num_pix[irr],sizeof(flouble));
    }
  }

#pragma omp parallel default(none)			\
  shared(oi,z0_node,zf_node,dx)
  {
    int ir;

#pragma omp for schedule(dynamic,DYNAMIC_SIZE)
    for(ir=0;ir<oi->nr;ir++) {
      double r0=oi->r0_arr[ir];
      double rf=oi->rf_arr[ir];
      double rm=0.5*(r0+rf);
      long ipx,npix=12*oi->nside_arr[ir]*oi->nside_arr[ir];
      flouble dr_trans=rm*sqrt(4*M_PI/npix);
      double dz_additional=dr_trans+dx;
      double z0_here=z0_node-dz_additional;
      double zf_here=zf_node+dz_additional;

      oi->num_pix[ir]=0;
      for(ipx=0;ipx<npix;ipx++) {
	double zpix0,zpixf,pos[3];
	//TODO: is this better than pix2ang?
	pix2vec_ring(oi->nside_arr[ir],ipx,pos);
	zpix0=r0*pos[2];
	zpixf=rf*pos[2];
	if(((zpix0>=z0_here) && (zpix0<=zf_here)) ||
	   ((zpixf>=z0_here) && (zpixf<=zf_here))) {
	  oi->list_ipix[ir][oi->num_pix[ir]]=ipx;
	  oi->num_pix[ir]++;
	}
      }
    } //end omp for
  } //end omp parallel

#ifdef _DEBUG
  for(irr=0;irr<oi->nr;irr++) {
    double r0=oi->r0_arr[irr];
    double rf=oi->rf_arr[irr];
    double rm=0.5*(r0+rf);
    int nside_here=oi->nside_arr[irr];
    fprintf(par->f_dbg,
	    "  Shell %d, r=%lf, nside=%d, angular resolution %lf Mpc/h, cell size %lf, %d pixels\n",
	    irr,rm,nside_here,rm*sqrt(4*M_PI/(12*nside_here*nside_here)),dr,oi->num_pix[irr]);
  }
#endif //_DEBUG

  return oi;
}

void free_onion_info(OnionInfo *oi)
{
  if(oi->nr>0) {
    int ir;
    free(oi->r0_arr);
    free(oi->rf_arr);
    free(oi->nside_arr);
    free(oi->num_pix);
    for(ir=0;ir<oi->nr;ir++) {
      if(oi->num_pix[ir]>0) {
	free(oi->list_ipix[ir]);
	free(oi->maps[ir]);
      }
    }
    free(oi->list_ipix);
    free(oi->maps);
  }
  free(oi);
}
