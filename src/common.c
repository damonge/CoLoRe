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

static inline flouble get_res(int nside)
{
#if PIXTYPE == PT_CEA
  return acos(1-2.0/nside);
#elif PIXTYPE == PT_CAR
  return M_PI/nside;
#else
  return acos(1-2.0/nside);
#endif //PIXTYPE
  //  return sqrt(2*M_PI)/nside;
}

static inline flouble get_lat_index(int nside,flouble cth)
{
#if PIXTYPE == PT_CEA
  return 0.5*(cth+1)*nside;
#elif PIXTYPE == PT_CAR
  return (1-acos(CLAMP(cth,-1,1))/M_PI)*nside;
#else
  return 0.5*(cth+1)*nside;
#endif //PIXTYPE
}

OnionInfo *alloc_onion_empty(ParamCoLoRe *par,int nside_base)
{
  int ir;
  OnionInfo *oi=my_malloc(sizeof(OnionInfo));
  double dx=par->l_box/par->n_grid;
  double dr=FAC_CART2SPH_VOL*dx;

  oi->nr=(int)(par->r_max/dr+1);
  dr=par->r_max/oi->nr; //Redefine radial interval
  oi->r0_arr=my_malloc(oi->nr*sizeof(flouble));
  oi->rf_arr=my_malloc(oi->nr*sizeof(flouble));
  oi->nside_arr=my_malloc(oi->nr*sizeof(int));
  oi->iphi0_arr=my_malloc(oi->nr*sizeof(int));
  oi->iphif_arr=my_malloc(oi->nr*sizeof(int));
  oi->icth0_arr=my_malloc(oi->nr*sizeof(int));
  oi->icthf_arr=my_malloc(oi->nr*sizeof(int));
  oi->num_pix=my_malloc(oi->nr*sizeof(int));

  for(ir=0;ir<oi->nr;ir++) {
    int nside_here=nside_base;
    flouble rm=(ir+0.5)*dr;
    flouble dr_trans=rm*get_res(nside_here);
    oi->r0_arr[ir]=ir*dr;
    oi->rf_arr[ir]=(ir+1)*dr;

    while(dr_trans>FAC_CART2SPH_VOL*dx) {
      nside_here*=2;
      dr_trans=rm*get_res(nside_here);
    }
    oi->nside_arr[ir]=nside_here;
  }

  return oi;
}

OnionInfo *alloc_onion_info_slices(ParamCoLoRe *par)
{
  int ir;
  double dx=par->l_box/par->n_grid;
  double z0_node=dx*par->iz0_here-par->pos_obs[2];
  double zf_node=dx*(par->iz0_here+par->nz_here)-par->pos_obs[2];
  OnionInfo *oi=alloc_onion_empty(par,par->nside_base);

  for(ir=0;ir<oi->nr;ir++) {
    int icth0,icthf;
    double cth0,cthf;
    double r0=oi->r0_arr[ir];
    double rf=oi->rf_arr[ir];
    double rm=0.5*(r0+rf);
    long npix=2*oi->nside_arr[ir]*oi->nside_arr[ir];
    flouble dr_trans=rm*get_res(oi->nside_arr[ir]);
    double dz_additional=dr_trans+dx;
    double z0_here=z0_node-dz_additional;
    double zf_here=zf_node+dz_additional;
    
    //Select phi bounds
    //Pick the whole ring in this case
    oi->iphi0_arr[ir]=0;
    oi->iphif_arr[ir]=2*oi->nside_arr[ir]-1;
    
    //Select cth bounds
    if(r0>0) {
      cth0=MIN(z0_here/rf,z0_here/r0);
      cthf=MAX(zf_here/rf,zf_here/r0);
    }
    else {
      cth0=z0_here/rf;
      cthf=zf_here/rf;
    }
    icth0=(int)(get_lat_index(oi->nside_arr[ir],cth0)-1);
    icthf=(int)(get_lat_index(oi->nside_arr[ir],cthf)+1);
    icth0=CLAMP(icth0,0,oi->nside_arr[ir]-1);
    icthf=CLAMP(icthf,0,oi->nside_arr[ir]-1);
    oi->icth0_arr[ir]=icth0;
    oi->icthf_arr[ir]=icthf;

    if(((oi->icth0_arr[ir]==oi->nside_arr[ir]-1) &&
	(oi->icthf_arr[ir]==oi->nside_arr[ir]-1) &&
	(z0_here>rf)) ||
       ((oi->icth0_arr[ir]==0) && (oi->icthf_arr[ir]==0) && (zf_here<-rf)))
      oi->num_pix[ir]=0;
    else
      oi->num_pix[ir]=2*oi->nside_arr[ir]*(icthf-icth0+1);

#ifdef _DEBUG
    int nside_here=oi->nside_arr[ir];
    fprintf(par->f_dbg,
	    "  Shell %d, r=%lf, nside=%d, angular resolution %lf Mpc/h, cell size %lf, %d pixels\n",
	    ir,rm,nside_here,rm*get_res(nside_here),dx,oi->num_pix[ir]);
#endif //_DEBUG
  }

  return oi;
}

OnionInfo **alloc_onion_info_beams(ParamCoLoRe *par)
{
  int i_base,ir,i_base_here,icth;
  OnionInfo **oi;
  OnionInfo *oi_dum=alloc_onion_empty(par,par->nside_base);
  int nside_base=oi_dum->nside_arr[0];
  int npix_base=2*nside_base*nside_base;
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

  i_base=0; i_base_here=0;
  for(icth=0;icth<nside_base;icth++) {
    int iphi;
    for(iphi=0;iphi<2*nside_base;iphi++) {
      if(i_base%NNodes==NodeThis) {
	for(ir=0;ir<oi[i_base_here]->nr;ir++) {
	  int nside_ratio=oi[i_base_here]->nside_arr[ir]/nside_base;
	  oi[i_base_here]->iphi0_arr[ir]=iphi*nside_ratio;
	  oi[i_base_here]->iphif_arr[ir]=(iphi+1)*nside_ratio-1;
	  oi[i_base_here]->icth0_arr[ir]=icth*nside_ratio;
	  oi[i_base_here]->icthf_arr[ir]=(icth+1)*nside_ratio-1;
	  oi[i_base_here]->num_pix[ir]=nside_ratio*nside_ratio;
	}
	i_base_here++;
      }
      i_base++;
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
    free(oi->iphi0_arr);
    free(oi->iphif_arr);
    free(oi->icth0_arr);
    free(oi->icthf_arr);
    free(oi->num_pix);
  }
  free(oi);
}

void free_slices(ParamCoLoRe *par)
{
  int ii;
  for(ii=0;ii<par->oi_slices->nr;ii++) {
    free(par->dens_slices[ii]);
    free(par->vrad_slices[ii]);
    if(par->do_lensing) {
      free(par->p_xx_slices[ii]);
      free(par->p_xy_slices[ii]);
      free(par->p_yy_slices[ii]);
    }
  }
  free(par->dens_slices);
  free(par->vrad_slices);
  if(par->do_lensing) {
    free(par->p_xx_slices);
    free(par->p_xy_slices);
    free(par->p_yy_slices);
  }
}

void alloc_slices(ParamCoLoRe *par)
{
  int ii;

  par->dens_slices=my_malloc(par->oi_slices->nr*sizeof(flouble *)); 
  par->vrad_slices=my_malloc(par->oi_slices->nr*sizeof(flouble *));
  if(par->do_lensing) {
    par->p_xx_slices=my_malloc(par->oi_slices->nr*sizeof(flouble *));
    par->p_xy_slices=my_malloc(par->oi_slices->nr*sizeof(flouble *));
    par->p_yy_slices=my_malloc(par->oi_slices->nr*sizeof(flouble *));
  }
  for(ii=0;ii<par->oi_slices->nr;ii++) {
    par->dens_slices[ii]=my_calloc(par->oi_slices->num_pix[ii],sizeof(flouble));
    par->vrad_slices[ii]=my_calloc(par->oi_slices->num_pix[ii],sizeof(flouble));
    if(par->do_lensing) {
      par->p_xx_slices[ii]=my_calloc(par->oi_slices->num_pix[ii],sizeof(flouble));
      par->p_xy_slices[ii]=my_calloc(par->oi_slices->num_pix[ii],sizeof(flouble));
      par->p_yy_slices[ii]=my_calloc(par->oi_slices->num_pix[ii],sizeof(flouble));
    }
  }
}

void free_beams(ParamCoLoRe *par)
{
  int ib;

  for(ib=0;ib<par->n_beams_here;ib++) {
    int ii;
    for(ii=0;ii<par->oi_slices->nr;ii++) {
      free(par->dens_beams[ib][ii]);
      free(par->vrad_beams[ib][ii]);
      if(par->do_lensing) {
	free(par->p_xx_beams[ib][ii]);
	free(par->p_xy_beams[ib][ii]);
	free(par->p_yy_beams[ib][ii]);
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
    free(par->nsrc_beams[ib]);
  }
  free(par->dens_beams);
  free(par->vrad_beams);
  if(par->do_lensing) {
    free(par->p_xx_beams);
    free(par->p_xy_beams);
    free(par->p_yy_beams);
  }
  free(par->nsrc_beams);
}

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
    par->nsrc_beams[ib]=my_malloc(par->oi_beams[ib]->nr*sizeof(int *));
    for(ii=0;ii<par->oi_slices->nr;ii++) {
      par->dens_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
      par->vrad_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
      if(par->do_lensing) {
	par->p_xx_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
	par->p_xy_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
	par->p_yy_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(flouble));
      }
      par->nsrc_beams[ib][ii]=my_calloc(par->oi_beams[ib]->num_pix[ii],sizeof(int));
    }
  }
}

void free_hp_shell(HealpixShells *shell)
{
  free(shell->listpix);
  free(shell->r0);
  free(shell->rf);
  free(shell->data);
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
