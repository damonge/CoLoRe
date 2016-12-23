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

flouble get_res(int nside)
{
#if PIXTYPE == PT_CEA
  return sqrt(2*M_PI)/nside;
  //  return acos(1-2.0/nside);
#elif PIXTYPE == PT_CAR
  return M_PI/nside;
#elif PIXTYPE == PT_HPX
  return sqrt(M_PI/3)/nside;
#else
  return sqrt(2*M_PI)/nside;
#endif //PIXTYPE
}

long get_npix(int nside)
{
#if PIXTYPE == PT_CEA
  return 2*nside*((long)nside);
#elif PIXTYPE == PT_CAR
  return 2*nside*((long)nside);
#elif PIXTYPE == PT_HPX
  return 12*nside*((long)nside);
#else
  return 2*nside*((long)nside);
#endif //PIXTYPE
}

void get_vec(int ipix_nest,int iphi_0,int icth_0,int nside,int nside_ratio,double *u)
{

#if PIXTYPE == PT_HPX
  pix2vec_nest(nside,iphi_0+ipix_nest,u);
#else

  double cth,sth,phi;
  int icth=ipix_nest/nside_ratio;
  int iphi=ipix_nest-icth*nside_ratio;

#if PIXTYPE == PT_CEA
  cth=(icth+icth_0+0.5)*2./nside-1;
#else PIXTYPE == PT_CAR
  cth=cos(M_PI-(icth+icth_0+0.5)*M_PI/nside);
#endif //PIXTYPE

  phi=(iphi+iphi_0+0.5)*M_PI/nside;
  sth=sqrt((1+cth)*(1-cth));
  u[0]=sth*cos(phi);
  u[1]=sth*sin(phi);
  u[2]=cth;

#endif //PIXTYPE
}

#define NS_RANDOM_EXTRA_HPX 4
void get_random_angles(gsl_rng *rng,int ipix_nest,int iphi_0,int icth_0,int nside,int nside_ratio,
		       double *th,double *phi)
{
#if PIXTYPE == PT_HPX
  int i_extra=(int)(NS_RANDOM_EXTRA_HPX*NS_RANDOM_EXTRA_HPX*rng_01(rng));
  pix2ang_nest(nside*NS_RANDOM_EXTRA_HPX,
	       (iphi_0+ipix_nest)*NS_RANDOM_EXTRA_HPX*NS_RANDOM_EXTRA_HPX+i_extra,th,phi);
  (*phi)+=(rng_01(rng)-0.5)*0.57/(nside*NS_RANDOM_EXTRA_HPX);
  (*th)+=(rng_01(rng)-0.5)*0.57/(nside*NS_RANDOM_EXTRA_HPX);
#else 

  int icth=ipix_nest/nside_ratio;
  int iphi=ipix_nest-icth*nside_ratio;

#if PIXTYPE == PT_CEA
  *th=acos(-1+(icth_0+icth+rng_01(rng))*2./nside);
#elif PIXTYPE == PT_CAR
  *th=M_PI-(icth_0+icth+rng_01(rng))*M_PI/nside;
#endif //PIXTYPE

  *phi=M_PI*(iphi_0+iphi+rng_01(rng))/nside;

#endif //PIXTYPE
}

OnionInfo *alloc_onion_empty(ParamCoLoRe *par,int nside_base)
{
  int ir;
  OnionInfo *oi=my_malloc(sizeof(OnionInfo));
  double dx=par->l_box/par->n_grid;
  double dr=FAC_CART2SPH_PAR*dx;

  oi->nr=(int)(par->r_max/dr+1);
  dr=par->r_max/oi->nr; //Redefine radial interval
  oi->r0_arr=my_malloc(oi->nr*sizeof(flouble));
  oi->rf_arr=my_malloc(oi->nr*sizeof(flouble));
  oi->nside_arr=my_malloc(oi->nr*sizeof(int));
  oi->nside_ratio_arr=my_malloc(oi->nr*sizeof(int));
  oi->iphi0_arr=my_malloc(oi->nr*sizeof(int));
  oi->icth0_arr=my_malloc(oi->nr*sizeof(int));
  oi->num_pix=my_malloc(oi->nr*sizeof(int));

  for(ir=0;ir<oi->nr;ir++) {
    int nside_here=nside_base;
    flouble rm=(ir+0.5)*dr;
    flouble dr_trans=rm*get_res(nside_here);
    oi->r0_arr[ir]=ir*dr;
    oi->rf_arr[ir]=(ir+1)*dr;

    while(dr_trans>FAC_CART2SPH_PERP*dx) {
      nside_here*=2;
      dr_trans=rm*get_res(nside_here);
    }
    oi->nside_arr[ir]=nside_here;
  }

  return oi;
}

OnionInfo **alloc_onion_info_beams(ParamCoLoRe *par)
{
  int i_base,ir,i_base_here;
  OnionInfo **oi;
  OnionInfo *oi_dum=alloc_onion_empty(par,par->nside_base);
  int nside_base=oi_dum->nside_arr[0];
  int npix_base=get_npix(nside_base);
  int nbase_per_node=npix_base/NNodes;
  int nbase_extra=npix_base%NNodes;
  int nbase_here=nbase_per_node;
  if(NodeThis<nbase_extra)
    nbase_here++;
  free_onion_info(oi_dum);

  par->n_beams_here=nbase_here;

  oi=my_malloc(nbase_here*sizeof(OnionInfo *));
  for(i_base=0;i_base<nbase_here;i_base++)
    oi[i_base]=alloc_onion_empty(par,nside_base);

  i_base_here=0;
  for(i_base=0;i_base<npix_base;i_base++) {
#if PIXTYPE != PT_HPX
    int icth=i_base/(2*nside_base);
    int iphi=i_base-icth*2*nside_base;
#endif //PIXTYPE
    if(i_base%NNodes==NodeThis) {
      for(ir=0;ir<oi[i_base_here]->nr;ir++) {
	int nside_ratio=oi[i_base_here]->nside_arr[ir]/nside_base;
	  oi[i_base_here]->nside_ratio_arr[ir]=nside_ratio;
#if PIXTYPE == PT_CEA || PIXTYPE == PT_CAR
	oi[i_base_here]->iphi0_arr[ir]=iphi*nside_ratio;
	oi[i_base_here]->icth0_arr[ir]=icth*nside_ratio;
#elif PIXTYPE == PT_HPX
	oi[i_base_here]->iphi0_arr[ir]=i_base*nside_ratio*nside_ratio;
	oi[i_base_here]->icth0_arr[ir]=-1;
#endif //PIXTYPE
	oi[i_base_here]->num_pix[ir]=nside_ratio*nside_ratio;
      }
      i_base_here++;
    }
  }

#ifdef _DEBUG
  double dx=par->l_box/par->n_grid;
  fprintf(par->f_dbg,"Beams: %d\n",par->n_beams_here);
  for(ir=0;ir<oi[0]->nr;ir++) {
    int ib;
    double r0=oi[0]->r0_arr[ir];
    double rf=oi[0]->rf_arr[ir];
    double rm=0.5*(r0+rf);
    int nside_here=oi[0]->nside_arr[ir];
    fprintf(par->f_dbg,
	    "  Shell %d, r=%lf, nside=%d, angular resolution %lf Mpc/h, cell size %lf",
	    ir,rm,nside_here,rm*get_res(nside_here),dx);
    fprintf(par->f_dbg,"[ ");
    for(ib=0;ib<par->n_beams_here;ib++)
      fprintf(par->f_dbg,"%d ",oi[ib]->num_pix[ir]);
    fprintf(par->f_dbg,"]\n");
  }
#endif //_DEBUG

  return oi;
}

void free_onion_info(OnionInfo *oi)
{
  if(oi->nr>0) {
    free(oi->r0_arr);
    free(oi->rf_arr);
    free(oi->nside_arr);
    free(oi->nside_ratio_arr);
    free(oi->iphi0_arr);
    free(oi->icth0_arr);
    free(oi->num_pix);
  }
  free(oi);
}

void free_beams(ParamCoLoRe *par)
{
  int ib;

  for(ib=0;ib<par->n_beams_here;ib++) {
    int ii;
    for(ii=0;ii<par->oi_beams[ib]->nr;ii++) {
      free(par->dens_beams[ib][ii]);
      free(par->vrad_beams[ib][ii]);
      if(par->do_lensing) {
	free(par->p_xx_beams[ib][ii]);
	free(par->p_xy_beams[ib][ii]);
	free(par->p_yy_beams[ib][ii]);
      }
      if(par->do_isw) {
	free(par->pdot_beams[ib][ii]);
      }
      free(par->nsrc_beams[ib][ii]);
    }
    free(par->dens_beams[ib]);
    free(par->vrad_beams[ib]);
    if(par->do_lensing) {
      free(par->p_xx_beams[ib]);
      free(par->p_xy_beams[ib]);
      free(par->p_yy_beams[ib]);
    }
    if(par->do_isw) {
      free(par->pdot_beams[ib]);
    }
    free(par->nsrc_beams[ib]);
  }
  free(par->dens_beams);
  free(par->vrad_beams);
  if(par->do_lensing) {
    free(par->p_xx_beams);
    free(par->p_xy_beams);
    free(par->p_yy_beams);
  }
  if(par->do_isw) {
    free(par->pdot_beams);
  }
  free(par->nsrc_beams);
}

unsigned long long get_max_memory(ParamCoLoRe *par)
{
  int ib;
  unsigned long long total_GB=0;
  unsigned long long total_GB_pix=0;
  unsigned long long total_GB_lpt=0;

  if(par->need_onions) {
    for(ib=0;ib<par->n_beams_here;ib++) {
      int ii;
      for(ii=0;ii<par->oi_beams[ib]->nr;ii++) {
	total_GB_pix+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
	total_GB_pix+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
	if(par->do_lensing) {
	  total_GB_pix+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
	  total_GB_pix+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
	  total_GB_pix+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
	}
	if(par->do_isw)
	  total_GB_pix+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
	total_GB_pix+=par->oi_beams[ib]->num_pix[ii]*sizeof(int);
      }
    }
  }

  if(par->dens_type==DENS_TYPE_1LPT) {
    total_GB_lpt=(unsigned long long)(3*(1+par->lpt_buffer_fraction)*par->nz_here*
				      ((lint)(par->n_grid*par->n_grid))*6*sizeof(flouble));
  }
  else if(par->dens_type==DENS_TYPE_2LPT) {
    total_GB_lpt=0;
    total_GB_lpt=(unsigned long long)(8*(1+par->lpt_buffer_fraction)*par->nz_here*
				      ((lint)(par->n_grid*par->n_grid))*6*sizeof(flouble));
  }

  total_GB=MAX(total_GB_lpt,total_GB_pix);

  void *ptest=my_malloc(total_GB);
  free(ptest);
  return total_GB;
}

void alloc_beams(ParamCoLoRe *par)
{
  int ib;
#ifdef _DEBUG
  unsigned long long total_GB=0;
  double total_GB_d=0;
#endif //_DEBUG

  par->dens_beams=my_malloc(par->n_beams_here*sizeof(flouble **));
  par->vrad_beams=my_malloc(par->n_beams_here*sizeof(flouble **));
  if(par->do_lensing) {
    par->p_xx_beams=my_malloc(par->n_beams_here*sizeof(flouble **));
    par->p_xy_beams=my_malloc(par->n_beams_here*sizeof(flouble **));
    par->p_yy_beams=my_malloc(par->n_beams_here*sizeof(flouble **));
  }
  if(par->do_isw) {
    par->pdot_beams=my_malloc(par->n_beams_here*sizeof(flouble **));
  }
  par->nsrc_beams=my_malloc(par->n_beams_here*sizeof(int **));
  for(ib=0;ib<par->n_beams_here;ib++) {
    int ii;
    par->dens_beams[ib]=my_malloc(par->oi_beams[ib]->nr*sizeof(flouble *));
    par->vrad_beams[ib]=my_malloc(par->oi_beams[ib]->nr*sizeof(flouble *));
    if(par->do_lensing) {
      par->p_xx_beams[ib]=my_malloc(par->oi_beams[ib]->nr*sizeof(flouble *));
      par->p_xy_beams[ib]=my_malloc(par->oi_beams[ib]->nr*sizeof(flouble *));
      par->p_yy_beams[ib]=my_malloc(par->oi_beams[ib]->nr*sizeof(flouble *));
    }
    if(par->do_isw) {
      par->pdot_beams[ib]=my_malloc(par->oi_beams[ib]->nr*sizeof(flouble *));
    }
    par->nsrc_beams[ib]=my_malloc(par->oi_beams[ib]->nr*sizeof(int *));
    for(ii=0;ii<par->oi_beams[ib]->nr;ii++) {
      par->dens_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
      par->vrad_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
#ifdef _DEBUG
      total_GB+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
      total_GB+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
#endif //_DEBUG
      if(par->do_lensing) {
	par->p_xx_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
	par->p_xy_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
	par->p_yy_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
#ifdef _DEBUG
	total_GB+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
	total_GB+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
	total_GB+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
#endif //_DEBUG
      }
      if(par->do_isw) {
	par->pdot_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
#ifdef _DEBUG
	total_GB+=par->oi_beams[ib]->num_pix[ii]*sizeof(flouble);
#endif //_DEBUG
      }
      par->nsrc_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(int));
#ifdef _DEBUG
      total_GB+=par->oi_beams[ib]->num_pix[ii]*sizeof(int);
#endif //_DEBUG
    }
  }

#ifdef _DEBUG
  total_GB_d=total_GB/pow(1024.,3);
  printf(" Node %d: Have allocated %.3lf GB for beams\n",NodeThis,total_GB_d);
#endif //_DEBUG
}

void free_hp_shell(HealpixShells *shell)
{
  free(shell->listpix);
  free(shell->r0);
  free(shell->rf);
  free(shell->data);
  free(shell->nadd);
}

HealpixShells *new_hp_shell(int nside,int nr)
{
  HealpixShells *shell=my_malloc(sizeof(HealpixShells));
  shell->nside=nside;

  //Figure out radial shells
  shell->nr=nr;
  shell->r0=my_malloc(shell->nr*sizeof(flouble));
  shell->rf=my_malloc(shell->nr*sizeof(flouble));

  return shell;
}
