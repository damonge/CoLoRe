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
  print_info(" MPIThreads = %d\n",MPIThreadsOK);
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

/*
#define NS_RANDOM_EXTRA_HPX 4
void get_random_angles(gsl_rng *rng,int ipix_nest,int ipix0,int nside,double *th,double *phi)
{
  long n_extra=NS_RANDOM_EXTRA_HPX;
  while(n_extra*nside>NSIDE_MAX_HPX)
    n_extra/=2;
  if(n_extra<=0) n_extra=1; //This should never happen
  
  int i_extra=(int)(n_extra*n_extra*rng_01(rng));
  pix2ang_nest(nside*n_extra,
	       (ipix0+ipix_nest)*n_extra*n_extra+i_extra,th,phi);
  (*phi)+=(rng_01(rng)-0.5)*0.57/(nside*n_extra);
  (*th)+=(rng_01(rng)-0.5)*0.57/(nside*n_extra);
}

static OnionInfo *alloc_onion_empty(ParamCoLoRe *par,int nside_base)
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
  oi->ipix0_arr=my_malloc(oi->nr*sizeof(int));
  oi->num_pix=my_malloc(oi->nr*sizeof(int));

  for(ir=0;ir<oi->nr;ir++) {
    int nside_here=nside_base;
    flouble rm=(ir+0.5)*dr;
    flouble dr_trans=rm*sqrt(M_PI/3)/nside_here;
    oi->r0_arr[ir]=ir*dr;
    oi->rf_arr[ir]=(ir+1)*dr;

    while((dr_trans>FAC_CART2SPH_PERP*dx) && (nside_here<NSIDE_MAX_HPX)) {
      nside_here*=2;
      dr_trans=rm*sqrt(M_PI/3)/nside_here;
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
  int npix_base=he_nside2npix(nside_base);
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
    if(i_base%NNodes==NodeThis) {
      for(ir=0;ir<oi[i_base_here]->nr;ir++) {
	int nside_ratio=oi[i_base_here]->nside_arr[ir]/nside_base;
	oi[i_base_here]->nside_ratio_arr[ir]=nside_ratio;
	oi[i_base_here]->ipix0_arr[ir]=i_base*nside_ratio*nside_ratio;
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
	    ir,rm,nside_here,rm*sqrt(M_PI/3)/nside_here,dx);
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
    free(oi->ipix0_arr);
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
*/

unsigned long long get_max_memory(ParamCoLoRe *par)
{
  unsigned long long total_GB=0;
  unsigned long long total_GB_gau=0;
  unsigned long long total_GB_pix=0;
  unsigned long long total_GB_lpt=0;

  total_GB_gau=(2*(par->nz_here+1)*((long)(par->n_grid*(par->n_grid/2+1))))*sizeof(dftw_complex);

  /*
  if(par->need_onions) {
    int ib;
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
  */

  if(par->dens_type==DENS_TYPE_1LPT) {
    total_GB_lpt=(unsigned long long)(3*(1+par->lpt_buffer_fraction)*par->nz_here*
				      ((long)((par->n_grid/2+1)*par->n_grid))*sizeof(dftw_complex));
  }
  else if(par->dens_type==DENS_TYPE_2LPT) {
    total_GB_lpt=0;
    total_GB_lpt=(unsigned long long)(8*(1+par->lpt_buffer_fraction)*par->nz_here*
				      ((long)((par->n_grid/2+1)*par->n_grid))*sizeof(dftw_complex));
  }

  total_GB=total_GB_gau+MAX(total_GB_lpt,total_GB_pix);

#ifdef _DEBUG
  printf("Node %d will allocate %.3lf GB (Gaussian)",NodeThis,(double)(total_GB_gau/pow(1024.,3)));
  if((par->dens_type==DENS_TYPE_1LPT) || (par->dens_type==DENS_TYPE_2LPT))
    printf(", %.3lf GB (%dLPT)",(double)(total_GB_lpt/pow(1024.,3)),par->dens_type);
  //  if(par->need_onions)
  //    printf(", %.3lf GB (pixels)",(double)(total_GB_pix/pow(1024.,3)));
  printf("\n");
#endif //_DEBUG

  void *ptest=my_malloc(total_GB);
  free(ptest);
  return total_GB;
}

void get_radial_params(double rmax,int ngrid,int *nr,double *dr)
{
  *nr=NSAMP_RAD*ngrid/2;
  *dr=rmax/(*nr);
}

/*
void alloc_beams(ParamCoLoRe *par)
{
  int ib;

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
      if(par->do_lensing) {
	par->p_xx_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
	par->p_xy_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
	par->p_yy_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
      }
      if(par->do_isw) {
	par->pdot_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
      }
      par->nsrc_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(int));
    }
  }
}
*/

CatalogCartesian *catalog_cartesian_alloc(int nsrcs)
{
  CatalogCartesian *cat=my_malloc(sizeof(CatalogCartesian));
  
  if(nsrcs>0) {
    cat->nsrc=nsrcs;
    cat->pos=my_malloc(NPOS_CC*nsrcs*sizeof(float));
    cat->ipix=my_malloc(nsrcs*sizeof(int));
  }
  else {
    cat->nsrc=0;
    cat->pos=NULL;
    cat->ipix=NULL;
  }

  return cat;
}

void catalog_cartesian_free(CatalogCartesian *cat)
{
  if(cat->nsrc>0) {
    free(cat->pos);
    free(cat->ipix);
  }
  free(cat);
}

Catalog *catalog_alloc(int nsrcs,int has_skw,double rmax,int ng)
{
  Catalog *cat=my_malloc(sizeof(Catalog));
  
  if(nsrcs>0) {
    cat->nsrc=nsrcs;
    cat->srcs=my_malloc(nsrcs*sizeof(Src));
    cat->has_skw=has_skw;
    get_radial_params(rmax,ng,&(cat->nr),&(cat->dr));
    cat->rmax=rmax;
    cat->idr=1./cat->dr;
    if(has_skw) {
      cat->d_skw=my_calloc(cat->nsrc*cat->nr,sizeof(float));
      cat->v_skw=my_calloc(cat->nsrc*cat->nr,sizeof(float));
    }
  }
  else {
    cat->nsrc=0;
    cat->srcs=NULL;
    cat->has_skw=0;
  }

  return cat;
}

void catalog_free(Catalog *cat)
{
  if(cat->nsrc>0) {
    free(cat->srcs);
    if(cat->has_skw) {
      free(cat->d_skw);
      free(cat->v_skw);
    }
  }
  free(cat);
}

HealpixShells *hp_shell_alloc(int nside,int nside_base,int nr)
{
  if(nside>NSIDE_MAX_HPX)
    report_error(1,"Can't go beyond nside=%d\n",NSIDE_MAX_HPX);
  if(nside<nside_base)
    report_error(1,"Can't go below nside=%d\n",nside_base);

  int ib;
  long ip;
  HealpixShells *shell=my_malloc(sizeof(HealpixShells));

  //Figure out pixel angular positions
  int nbases=he_nside2npix(nside_base);
  int nside_ratio=nside/nside_base;
  long npix_perbeam=nside_ratio*nside_ratio;
  int nbeams_here=0;
  
  for(ib=NodeThis;ib<nbases;ib+=NNodes)
    nbeams_here++;
  shell->nside=nside;
  shell->num_pix=nside_ratio*nside_ratio*nbeams_here;
  shell->listpix=my_malloc(shell->num_pix*sizeof(long));
  shell->pos=my_malloc(3*shell->num_pix*sizeof(double));

  double *u=shell->pos;
  long ipix=0;
  for(ib=NodeThis;ib<nbases;ib+=NNodes) {
    for(ip=0;ip<npix_perbeam;ip++) {
      long id_nest=ib*npix_perbeam+ip;
      shell->listpix[ipix]=id_nest;
      pix2vec_nest(shell->nside,id_nest,u);
      u+=3;
      ipix++;
    }
  }

  //Figure out radial shells
  shell->nr=nr;
  shell->r0=my_malloc(shell->nr*sizeof(flouble));
  shell->rf=my_malloc(shell->nr*sizeof(flouble));

  //Zero all data and clear
  shell->data=my_calloc(shell->nr*shell->num_pix,sizeof(flouble));
  shell->nadd=my_calloc(shell->nr*shell->num_pix,sizeof(int));
  
  return shell;
}

void hp_shell_free(HealpixShells *shell)
{
  if(shell->listpix!=NULL)
    free(shell->listpix);
  if(shell->pos!=NULL)
    free(shell->pos);
  if(shell->r0!=NULL)
    free(shell->r0);
  if(shell->rf!=NULL)
    free(shell->rf);
  if(shell->data!=NULL)
    free(shell->data);
  if(shell->nadd!=NULL)
    free(shell->nadd);
  free(shell);
}
