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
*/

void get_radial_params(double rmax,int ngrid,int *nr,double *dr)
{
  *nr=NSAMP_RAD*ngrid/2;
  *dr=rmax/(*nr);
}

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

Catalog *catalog_alloc(int nsrcs,int has_lensing,int has_skw,
                       int skw_gauss,double rmax,int ng)
{
  Catalog *cat=my_malloc(sizeof(Catalog));

  if(nsrcs>0) {
    cat->nsrc=nsrcs;
    cat->srcs=my_malloc(nsrcs*sizeof(Src));
    cat->has_lensing=has_lensing;
    cat->has_skw=has_skw;
    cat->skw_gauss=skw_gauss;
    get_radial_params(rmax,ng,&(cat->nr),&(cat->dr));
    cat->rmax=rmax;
    cat->idr=1./cat->dr;
    if(has_skw) {
      if(skw_gauss) {
        cat->g_skw=my_calloc(cat->nsrc*cat->nr,sizeof(float));
      }
      else {
        cat->d_skw=my_calloc(cat->nsrc*cat->nr,sizeof(float));
      }
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
      if(cat->skw_gauss) {
        free(cat->g_skw);
      }
      else {
        free(cat->d_skw);
      }
      free(cat->v_skw);
    }
  }
  free(cat);
}

void hp_shell_adaptive_free(HealpixShellsAdaptive *shell)
{
  int ib,ir;
  free(shell->r);
  free(shell->nside);
  free(shell->num_pix_per_beam);
  for(ib=0;ib<shell->nbeams;ib++) {
    for(ir=0;ir<shell->nr;ir++)
      free(shell->data[ib][ir]);
    free(shell->data[ib]);
    free(shell->pos[ib]);
  }
  free(shell->data);
  free(shell->pos);
  free(shell);
}

static int *get_adaptive_nside(int nside_max, int nside_base,
                               int nr, flouble *r_arr, flouble dx, flouble dx_fraction)
{
  int ir;
  int *nsides=my_malloc(nr*sizeof(int));
  for(ir=0; ir<nr; ir++) {
    flouble r=r_arr[ir];
    int nside_here=nside_base;
    while(nside_here<nside_max) {
      double theta=sqrt(4*M_PI/he_nside2npix(nside_here));
      double dxt=theta*r;
      if(dxt<=dx*dx_fraction)
        break;
      nside_here*=2;
    }
    nsides[ir]=nside_here;
  }
  return nsides;
}
  
HealpixShellsAdaptive *hp_shell_adaptive_alloc(int nq, int nside_max, int nside_base,
                                               int nr, flouble *r_arr, flouble dx,
                                               flouble dx_fraction)
{
  if(nside_max>NSIDE_MAX_HPX)
    report_error(1,"Can't go beyond nside=%d\n",NSIDE_MAX_HPX);
  if(nside_max<nside_base)
    report_error(1,"Can't go below nside=%d\n",nside_base);
  
  HealpixShellsAdaptive *shell=my_malloc(sizeof(HealpixShellsAdaptive));

  int ib, ir;
  int nbases=he_nside2npix(nside_base);
  shell->nbeams=0;
  for(ib=NodeThis;ib<nbases;ib+=NNodes)
    shell->nbeams++;

  shell->nq=nq;
  shell->nr=nr;
  shell->nside=get_adaptive_nside(nside_max, nside_base, nr,
                                  r_arr, dx, dx_fraction);
  shell->num_pix_per_beam=my_malloc(nr*sizeof(long));
  shell->r=my_malloc(nr*sizeof(flouble));
  for(ir=0; ir<nr; ir++) {
    int nside_ratio;
    shell->r[ir]=r_arr[ir];
    nside_ratio=shell->nside[ir]/nside_base;
    shell->num_pix_per_beam[ir]=nside_ratio*nside_ratio;
  }

  long npix_hi=shell->num_pix_per_beam[nr-1];
  shell->data=my_malloc(shell->nbeams*sizeof(flouble **));
  shell->pos=my_malloc(shell->nbeams*sizeof(double **));
  for(ib=0;ib<shell->nbeams;ib++) {
    shell->data[ib]=my_malloc(nr*sizeof(flouble *));
    shell->pos[ib]=my_malloc(npix_hi*3*sizeof(double));
    for(ir=0; ir<nr; ir++)
      shell->data[ib][ir]=my_malloc(shell->nq*shell->num_pix_per_beam[ir]*sizeof(flouble));
  }

  for(ib=0;ib<shell->nbeams;ib++) {
    long ip, ip0=(ib*NNodes+NodeThis)*npix_hi;
    double *u=shell->pos[ib];
    for(ip=0;ip<npix_hi;ip++) {
      long id_nest=ip0+ip;
      pix2vec_nest(shell->nside[nr-1],id_nest,u);
      u+=3;
    }
  }
  return shell;
}
    
HealpixShells *hp_shell_alloc(int nq, int nside,int nside_base,int nr)
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
  shell->nq=nq;
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
  shell->data=my_calloc(shell->nq*shell->nr*shell->num_pix,sizeof(flouble));
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

unsigned long long get_max_memory(ParamCoLoRe *par,int just_test)
{
  unsigned long long total_GB=0;
  unsigned long long total_GB_gau=0;
  unsigned long long total_GB_lpt=0;
  int fac_gau=2;
  if(par->need_beaming) fac_gau=3;

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
    for(ipop=0;ipop<par->n_srcs;ipop++) {
      int nz,ii;
      long nsrc=0;
      double nztot=0;
      FILE *fi=fopen(par->fnameNzSrcs[ipop],"r");

      double *zarr,*nzarr;
      if(fi==NULL) error_open_file(par->fnameNzSrcs[ipop]);
      nz=linecount(fi); rewind(fi);
      zarr=my_malloc(nz*sizeof(double));
      nzarr=my_malloc(nz*sizeof(double));
      for(ii=0;ii<nz;ii++) {
	int stat=fscanf(fi,"%lf %lf",&(zarr[ii]),&(nzarr[ii]));
	if(stat!=2) error_read_line(par->fnameNzSrcs[ipop],ii+1);
	nzarr[ii]*=RTOD*RTOD;
      }
      for(ii=0;ii<nz-1;ii++) {
	if((zarr[ii]>=par->z_min) && (zarr[ii]<=par->z_max))
	  nztot+=nzarr[ii]*(zarr[ii+1]-zarr[ii]);
      }
      if((zarr[ii]>=par->z_min) && (zarr[ii]<=par->z_max))
	nztot+=nzarr[ii]*(zarr[ii+1]-zarr[ii]);
      nztot*=4*M_PI/NNodes;
      nsrc+=(long)(nztot);
      if(just_test)
	print_info(" Expect %ld type-%d sources\n",(long)(nztot*NNodes),ipop);
      free(zarr);
      free(nzarr);
      fclose(fi);

      long size_source=sizeof(Src)+NPOS_CC*sizeof(float)+sizeof(int);
      if(par->skw_srcs[ipop]) {
	int nr=NSAMP_RAD*par->n_grid/2;
	size_source+=2*nr*sizeof(float);
      }
      total_GB_srcs+=size_source*nsrc;
    }
  }

  unsigned long long total_GB_imap=0;
  if(par->do_imap) {
    int ipop;
    for(ipop=0;ipop<par->n_imap;ipop++) {
      int nr;
      unsigned long long size_map=he_nside2npix(par->nside_imap[ipop]);
      FILE *fnu=fopen(par->fnameNuImap[ipop],"r");
      if(fnu==NULL) error_open_file(par->fnameNuImap[ipop]);
      nr=linecount(fnu);
      fclose(fnu);
      total_GB_imap+=size_map*nr*(sizeof(flouble)+sizeof(int));
    }
  }

  unsigned long long total_GB_cstm=0;
  if(par->do_cstm) {
    int ipop;
    for(ipop=0;ipop<par->n_cstm;ipop++) {
      unsigned long long size_map=he_nside2npix(par->nside_cstm[ipop]);
      total_GB_cstm+=size_map*(sizeof(flouble)+sizeof(int));
    }
  }

  unsigned long long total_GB_kappa=0;
  if(par->do_kappa) {
    int ib;
    int nr=par->n_kappa;
    int nbases=he_nside2npix(par->nside_base);
    int nside_ratio=par->nside_kappa/par->nside_base;
    int npix_perbeam=nside_ratio*nside_ratio;
    unsigned long long size_map=0;
    for(ib=NodeThis;ib<nbases;ib+=NNodes)
      size_map+=npix_perbeam;
    total_GB_kappa+=size_map*nr*(sizeof(flouble)+sizeof(int));
  }

  unsigned long long total_GB_lensing=0;
  if(par->do_lensing) {
    int ib,ir;
    int nr=par->n_lensing;
    int nbases=he_nside2npix(par->nside_base);
    flouble *rs=compute_lensing_spacing(par);
    int *nsides=get_adaptive_nside(par->nside_lensing,
                                   par->nside_base, nr, rs,
                                   par->l_box/par->n_grid, 1.);
    int nbeams_here=0;
    for(ib=NodeThis;ib<nbases;ib+=NNodes)
      nbeams_here++;

    unsigned long long npix_total=0;
    for(ir=0;ir<nr;ir++) {
      int nside_ratio=nsides[ir]/par->nside_base;
      npix_total+=nside_ratio*nside_ratio;
    }
    npix_total*=nbeams_here;
    int nside_ratio_hi=nsides[nr-1]/par->nside_base;
    unsigned long long npix_hi=nside_ratio_hi*nside_ratio_hi*nbeams_here;
    total_GB_lensing+=npix_total*5*sizeof(flouble)+npix_hi*3*sizeof(double);
    free(rs);
    free(nsides);
  }

  unsigned long long total_GB_isw=0;
  if(par->do_isw) {
    int ib;
    int nr=par->n_isw;
    int nbases=he_nside2npix(par->nside_base);
    int nside_ratio=par->nside_isw/par->nside_base;
    int npix_perbeam=nside_ratio*nside_ratio;
    unsigned long long size_map=0;
    for(ib=NodeThis;ib<nbases;ib+=NNodes)
      size_map+=npix_perbeam;
    total_GB_isw+=size_map*nr*(sizeof(flouble)+sizeof(int));
  }
  total_GB=total_GB_gau+
    total_GB_lpt+
    total_GB_srcs+
    total_GB_imap+
    total_GB_cstm+
    total_GB_kappa+
    total_GB_lensing+
    total_GB_isw;

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
      if(par->do_imap)
	printf(", %.3lf GB (imap)",(double)(total_GB_imap/pow(1024.,3)));
      if(par->do_kappa)
	printf(", %.3lf MB (kappa)",(double)(total_GB_kappa/pow(1024.,2)));
      if(par->do_lensing)
	printf(", %.3lf MB (lensing)",(double)(total_GB_lensing/pow(1024.,2)));
      if(par->do_isw)
	printf(", %.3lf MB (isw)",(double)(total_GB_isw/pow(1024.,2)));
      if(par->do_cstm)
	printf(", %.3lf MB (custom)",(double)(total_GB_cstm/pow(1024.,2)));
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
