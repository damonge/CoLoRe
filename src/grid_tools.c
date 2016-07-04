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

#define N_SUBPART 10

static int get_inu(ParamCoLoRe *par,double nu,int inu_start,int ipop)
{
  int gotit=0;
  int inu0;
  if(inu_start<0)
    inu0=0;
  else if(inu_start>=par->n_nu[ipop])
    inu0=par->n_nu[ipop]-1;
  else
    inu0=inu_start;
  
  while(!gotit) {
    if((inu0==-1) || (inu0==par->n_nu[ipop]))
      gotit=1;
    else {
      if(nu<par->nu0_arr[ipop][inu0])
	inu0--;
      else {
	if(nu>=par->nuf_arr[ipop][inu0])
	  inu0++;
	else
	  gotit=1;
      }
    }
  }

  return inu0;
}

#ifdef _OLD_IM
static void get_IM_single(ParamCoLoRe *par,int ipop)
{
  //////
  // 1 - Applies lognormal transformation to density grid
  // 2 - Transforms radial velocity into redshift distortion
  // 3 - Substitutes the overdensity grid for the corresponding
  //     HI mass in each cell.
  int ii;
  double x_sub[N_SUBPART],y_sub[N_SUBPART],z_sub[N_SUBPART];
  long npix_ang=nside2npix(par->nside_im[ipop]);
  gsl_rng *rng=init_rng(par->seed_rng);
  double lcell=par->l_box/par->n_grid;
  for(ii=0;ii<N_SUBPART;ii++) {
    x_sub[ii]=lcell*(rng_01(rng)-0.5);
    y_sub[ii]=lcell*(rng_01(rng)-0.5);
    z_sub[ii]=lcell*(rng_01(rng)-0.5);
  }
  end_rng(rng);
  par->maps_IM[ipop]=my_calloc(par->n_nu[ipop]*npix_ang,sizeof(flouble));

  print_info(" %d-th IM species\n",ipop);
  if(NodeThis==0) timer(0);
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,ipop,npix_ang,x_sub,y_sub,z_sub)
#endif //_HAVE_OMP
  {
    lint iz;
    double dx=par->l_box/par->n_grid;
    int ngx=2*(par->n_grid/2+1);
    double r_max_shell=r_of_z(par,par->nu_rest[ipop]/par->nu0_arr[ipop][0]-1)+DR_RSD_ADDITIONAL;
    double r_min_shell=r_of_z(par,par->nu_rest[ipop]/par->nuf_arr[ipop][par->n_nu[ipop]-1]-1)-DR_RSD_ADDITIONAL;
    printf("%lf %lf\n",r_min_shell,r_max_shell);

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      int inu=0;
      lint indexz=iz*((lint)(ngx*par->n_grid));
      double z0=dx*(iz+par->iz0_here+0.5)-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint indexy=iy*ngx;
	double y0=dx*(iy+0.5)-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=ix+indexy+indexz;
	  double x0=dx*(ix+0.5)-par->pos_obs[0];
	  double rr=sqrt(x0*x0+y0*y0+z0*z0);
	  double redshift=z_of_r(par,rr);
	  if((rr<=r_max_shell) && (rr>=r_min_shell)) {
	    int isub;
	    double gfd=dgrowth_of_r(par,rr)*bias_of_z_im(par,redshift,ipop);
	    double dz_rsd=par->grid_rvel[index]*vgrowth_of_r(par,rr);
	    double dens_LN=exp(gfd*(par->grid_dens[index]-0.5*gfd*par->sigma2_gauss));
	    double mass_sub=temp_of_z_im(par,redshift,ipop)*dens_LN/N_SUBPART;
	    for(isub=0;isub<N_SUBPART;isub++) {
	      double x=x0+x_sub[isub];
	      double y=y0+y_sub[isub];
	      double z=z0+z_sub[isub];
	      double r=sqrt(x*x+y*y+z*z);
	      double zred=z_of_r(par,r)+dz_rsd;
	      double nu=par->nu_rest[ipop]/(1+zred);
	      inu=get_inu(par,nu,inu,ipop);
	      if((inu>=0)&&(inu<par->n_nu[ipop])) {
		long ipix;
		double pos[3]={x,y,z};

		vec2pix_ring(par->nside_im[ipop],pos,&ipix);
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
		par->maps_IM[ipop][ipix+npix_ang*inu]+=mass_sub;
	      }
	    }
	  }
	}
      }
    } //end pragma omp for
#ifdef _HAVE_OMP
#pragma omp barrier
#endif //_HAVE_OMP

    int inu;
#ifdef _HAVE_OMP
#pragma omp single 
#endif //_HAVE_OMP
    {
      print_info(" Normalizing to temperature\n");
    }
#ifdef _HAVE_OMP
#pragma omp for nowait
#endif //_HAVE_OMP
    for(inu=0;inu<par->n_nu[ipop];inu++) {
      long ipix;
      double pix_area=4*M_PI/npix_ang;
      double r0=r_of_z(par,par->nu_rest[ipop]/par->nuf_arr[ipop][inu]-1);
      double rf=r_of_z(par,par->nu_rest[ipop]/par->nu0_arr[ipop][inu]-1);
      double vol_ratio=dx*dx*dx/(pix_area*(rf*rf*rf-r0*r0*r0)*0.33333333333333);
      for(ipix=0;ipix<npix_ang;ipix++) {
	long index=ipix+npix_ang*inu;
	par->maps_IM[ipop][index]*=vol_ratio;
      }
    } //end pragma omp for
  } //end pragma omp parallel
  if(NodeThis==0) timer(2);

#ifdef _HAVE_MPI
#define REDUCE_BATCH 1073741824
  int remainder;
  long ipix0_here=0;
  long npix_all=par->n_nu[ipop]*nside2npix(par->nside_im[ipop]);
  while(ipix0_here+REDUCE_BATCH<=npix_all) {
    if(NodeThis==0) {
      MPI_Reduce(MPI_IN_PLACE,&(par->maps_IM[ipop][ipix0_here]),REDUCE_BATCH,
		 FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    }
    else {
      MPI_Reduce(&(par->maps_IM[ipop][ipix0_here]),NULL,REDUCE_BATCH,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    }
    ipix0_here+=REDUCE_BATCH;
  }
  remainder=(int)(npix_all-ipix0_here);
  if(remainder) {
    if(NodeThis==0) {
      MPI_Reduce(MPI_IN_PLACE,&(par->maps_IM[ipop][ipix0_here]),remainder,
		 FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    }
    else {
      MPI_Reduce(&(par->maps_IM[ipop][ipix0_here]),NULL,remainder,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    }
  }
#endif //_HAVE_MPI
}

#else //_OLD_IM

#ifdef _IM_3D22D
//TODO: there may be a faster way
static void find_pix_id(long ipix,long *list_ipix,int npix,int *id_inout)
{
  int found=0;
  int id_last=*id_inout;

  if(ipix>list_ipix[id_last]) {
    while((!found) && (id_last<npix)) {
      id_last++;
      if(ipix==list_ipix[id_last])
	found=1;
    }
    //    if(id_last==npix)
    //      report_error(1,"Pixel id not found!\n");
  }
  else if(ipix<list_ipix[id_last]) {
    while((!found) && (id_last>=0)) {
      id_last--;
      if(ipix==list_ipix[id_last])
	found=1;
    }
    //    if(id_last==-1)
    //      report_error(1,"Pixel id not found! %ld\n",ipix);
  }
  
  //  print_info("%d %d\n",*id_inout,id_last);

  *id_inout=id_last;
}

static void get_IM_single(ParamCoLoRe *par,int ipop)
{
  int ired;
  flouble *z0_arr=my_malloc(par->n_nu[ipop]*sizeof(flouble));
  flouble *zf_arr=my_malloc(par->n_nu[ipop]*sizeof(flouble));

  print_info(" %d-th IM species\n",ipop);
  if(NodeThis==0) timer(0);

  for(ired=0;ired<par->n_nu[ipop];ired++) {
    z0_arr[ired]=par->nu_rest[ipop]/par->nuf_arr[ipop][ired]-1;
    zf_arr[ired]=par->nu_rest[ipop]/par->nu0_arr[ipop][ired]-1;
  }
  //Collect pixels in this node.
  alloc_onion_info(par,&(par->oi_IM[ipop]),par->nside_im[ipop],
		   par->n_nu[ipop],z0_arr,zf_arr,1,0);
  free(z0_arr);
  free(zf_arr);

  int ii;
  long npix_ang=nside2npix(par->nside_im[ipop]);
  double x_sub[N_SUBPART],y_sub[N_SUBPART],z_sub[N_SUBPART];
  gsl_rng *rng=init_rng(par->seed_rng);
  double lcell=par->l_box/par->n_grid;
  for(ii=0;ii<N_SUBPART;ii++) {
    x_sub[ii]=lcell*(rng_01(rng)-0.5);
    y_sub[ii]=lcell*(rng_01(rng)-0.5);
    z_sub[ii]=lcell*(rng_01(rng)-0.5);
  }
  end_rng(rng);

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,ipop,npix_ang,x_sub,y_sub,z_sub)
#endif //_HAVE_OMP
  {
    lint iz;
    double dx=par->l_box/par->n_grid;
    int ngx=2*(par->n_grid/2+1);
    int *id_pix=my_calloc(par->n_nu[ipop],sizeof(int));
    double r_max_shell=r_of_z(par,par->nu_rest[ipop]/par->nu0_arr[ipop][0]-1)+DR_RSD_ADDITIONAL;
    double r_min_shell=r_of_z(par,par->nu_rest[ipop]/par->nuf_arr[ipop][par->n_nu[ipop]-1]-1)-DR_RSD_ADDITIONAL;
    printf("%lf %lf\n",r_min_shell,r_max_shell);

#ifdef _HAVE_OMP
#pragma omp for schedule(dynamic)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      int inu=0;
      lint indexz=iz*((lint)(ngx*par->n_grid));
      double z0=dx*(iz+par->iz0_here+0.5)-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint indexy=iy*ngx;
	double y0=dx*(iy+0.5)-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=ix+indexy+indexz;
	  double x0=dx*(ix+0.5)-par->pos_obs[0];
	  double rr=sqrt(x0*x0+y0*y0+z0*z0);
	  double redshift=z_of_r(par,rr);
	  if((rr<=r_max_shell) && (rr>=r_min_shell)) {
	    int isub;
	    double gfd=dgrowth_of_r(par,rr)*bias_of_z_im(par,redshift,ipop);
	    double dz_rsd=par->grid_rvel[index]*vgrowth_of_r(par,rr);
	    double dens_LN=exp(gfd*(par->grid_dens[index]-0.5*gfd*par->sigma2_gauss));
	    double mass_sub=temp_of_z_im(par,redshift,ipop)*dens_LN/N_SUBPART;
	    for(isub=0;isub<N_SUBPART;isub++) {
	      double x=x0+x_sub[isub];
	      double y=y0+y_sub[isub];
	      double z=z0+z_sub[isub];
	      double r=sqrt(x*x+y*y+z*z);
	      double zred=z_of_r(par,r)+dz_rsd;
	      double nu=par->nu_rest[ipop]/(1+zred);
	      inu=get_inu(par,nu,inu,ipop);
	      if((inu>=0)&&(inu<par->n_nu[ipop])) {
		long ipix;
		double pos[3]={x,y,z};
		long *lpx=(par->oi_IM[ipop]).list_ipix[inu];
		int npx=(par->oi_IM[ipop]).num_pix[inu];
		flouble *mps=(par->oi_IM[ipop]).maps[inu];

		vec2pix_ring(par->nside_im[ipop],pos,&ipix);
		find_pix_id(ipix,lpx,npx,&(id_pix[inu]));
		if((id_pix[inu]<0) || (id_pix[inu]>=npx)) {
		  double pospix[3];
		  double z_here;
		  double r0=r;
		  double dr=r_of_z(par,zred)-r0;
		  double dr_b=dz_rsd*ihub_of_r(par,r0);
		  pix2vec_ring(par->nside_im[ipop],ipix,pospix);
		  z_here=pospix[2]*(r0+dr);
		  report_error(1,"Wrong id, %d %ld. costh=%lf, r0=%lf, dr=%lf,%lf, z=%lf\n",
			       inu,ipix,pospix[2],r0,dr,dr_b,z_here);
		}
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
		mps[id_pix[inu]]+=mass_sub;
	      }
	    }
	  }
	}
      }
    } //end pragma omp for
    free(id_pix);
#ifdef _HAVE_OMP
#pragma omp barrier
#endif //_HAVE_OMP

    int inu;
#ifdef _HAVE_OMP
#pragma omp single 
#endif //_HAVE_OMP
    {
      print_info(" Normalizing to temperature\n");
    }
#ifdef _HAVE_OMP
#pragma omp for nowait
#endif //_HAVE_OMP
    for(inu=0;inu<par->n_nu[ipop];inu++) {
      long id;
      double pix_area=4*M_PI/npix_ang;
      double r0=r_of_z(par,par->nu_rest[ipop]/par->nuf_arr[ipop][inu]-1);
      double rf=r_of_z(par,par->nu_rest[ipop]/par->nu0_arr[ipop][inu]-1);
      double vol_ratio=dx*dx*dx/(pix_area*(rf*rf*rf-r0*r0*r0)*0.33333333333333);
      flouble *mps=(par->oi_IM[ipop]).maps[inu];
      for(id=0;id<(par->oi_IM[ipop]).num_pix[inu];id++) {
	mps[id]*=vol_ratio;
      }
    } //end pragma omp for
  } //end pragma omp parallel
  if(NodeThis==0) timer(2);
}

#else //_IM_3D22D

static void get_IM_single(ParamCoLoRe *par,int ipop)
{
  int ired;

  print_info(" %d-th IM species\n",ipop);
  if(NodeThis==0) timer(0);

  //Collect pixels in this node.
  flouble *z0_arr=my_malloc(par->n_nu[ipop]*sizeof(flouble));
  flouble *zf_arr=my_malloc(par->n_nu[ipop]*sizeof(flouble));
  for(ired=0;ired<par->n_nu[ipop];ired++) {
    z0_arr[ired]=par->nu_rest[ipop]/par->nuf_arr[ipop][ired]-1;
    zf_arr[ired]=par->nu_rest[ipop]/par->nu0_arr[ipop][ired]-1;
  }
  alloc_onion_info(par,&(par->oi_IM[ipop]),par->nside_im[ipop],
		   par->n_nu[ipop],z0_arr,zf_arr,1,0);
  free(z0_arr);
  free(zf_arr);
  
  //Compute redshift-space, light-cone, log-normal density field
  printf("Computing lognormal\n");
#pragma omp parallel default(none)		\
  shared(par,ipop)
  {
    lint iz;
    int ngx=2*(par->n_grid/2+1);
    double dx=par->l_box/par->n_grid;

#pragma omp for schedule(dynamic)
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint indexz=iz*((lint)(ngx*par->n_grid));
      double z0=dx*(iz+par->iz0_here+0.5)-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint indexy=iy*ngx;
	double y0=dx*(iy+0.5)-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=ix+indexy+indexz;
	  double x0=dx*(ix+0.5)-par->pos_obs[0];
	  double rr=sqrt(x0*x0+y0*y0+z0*z0);
	  double redshift=z_of_r(par,rr);
	  double gfd=dgrowth_of_r(par,rr)*bias_of_z_im(par,redshift,ipop);
	  //	  double dz_rsd=par->grid_rvel[index]*vgrowth_of_r(par,rr);
	  double dens_LN=exp(gfd*(par->grid_dens[index]-0.5*gfd*par->sigma2_gauss));
	  par->grid_vpot[index]=temp_of_z_im(par,redshift,ipop)*dens_LN;
	}
      }
    } //end omp for
  } //end omp parallel

  printf("Interpolating\n");
#pragma omp parallel default(none) \
  shared(par,ipop)
  {

    int inu;
    long nside=par->nside_im[ipop];
    double dx=par->l_box/par->n_grid;
    double i_dx=1./dx;
    int ngx=2*(par->n_grid/2+1);

#pragma omp for schedule(dynamic)
    for(inu=0;inu<par->n_nu[ipop];inu++) {
      int id;
      long *lpx=(par->oi_IM[ipop]).list_ipix[inu];
      int npx=(par->oi_IM[ipop]).num_pix[inu];
      flouble *mps=(par->oi_IM[ipop]).maps[inu];
      flouble r0=(par->oi_IM[ipop]).r0_arr[inu];
      flouble rf=(par->oi_IM[ipop]).rf_arr[inu];
      flouble rr=0.5*(r0+rf);
      printf("%d\n",inu);
      for(id=0;id<npx;id++) {
	int ax;
	lint i0[3],i1[3];
	double pos[3],u[3],v[3];
	int ipix=lpx[id];
	pix2vec_ring(nside,ipix,pos);
	for(ax=0;ax<3;ax++) {
	  pos[ax]=par->pos_obs[ax]+pos[ax]*rr-0.5*dx;
	  i0[ax]=(int)(pos[ax]*i_dx);
	  i1[ax]=i0[ax]+1;
	  u[ax]=pos[ax]*i_dx-i0[ax];
	  v[ax]=1-u[ax];
	}
	i0[2]-=par->iz0_here;
	i1[2]-=par->iz0_here;
	if((i0[2]>=0) && (i0[2]<par->nz_here) && (i0[1]>=0) && (i0[1]<par->n_grid) && (i0[0]>=0) && (i0[0]<par->n_grid))
	  mps[id]+=par->grid_vpot[i0[0]+ngx*(i0[1]+par->n_grid*i0[2])]*v[0]*v[1]*v[2];
	if((i0[2]>=0) && (i0[2]<par->nz_here) && (i0[1]>=0) && (i0[1]<par->n_grid) && (i1[0]>=0) && (i1[0]<par->n_grid))
	  mps[id]+=par->grid_vpot[i1[0]+ngx*(i0[1]+par->n_grid*i0[2])]*u[0]*v[1]*v[2];
	if((i0[2]>=0) && (i0[2]<par->nz_here) && (i1[1]>=0) && (i1[1]<par->n_grid) && (i0[0]>=0) && (i0[0]<par->n_grid))
	  mps[id]+=par->grid_vpot[i0[0]+ngx*(i1[1]+par->n_grid*i0[2])]*v[0]*u[1]*v[2];
	if((i0[2]>=0) && (i0[2]<par->nz_here) && (i1[1]>=0) && (i1[1]<par->n_grid) && (i1[0]>=0) && (i1[0]<par->n_grid))
	  mps[id]+=par->grid_vpot[i1[0]+ngx*(i1[1]+par->n_grid*i0[2])]*u[0]*u[1]*v[2];
	if((i1[2]>=0) && (i0[2]<par->nz_here) && (i0[1]>=0) && (i0[1]<par->n_grid) && (i0[0]>=0) && (i0[0]<par->n_grid))
	  mps[id]+=par->grid_vpot[i0[0]+ngx*(i0[1]+par->n_grid*i1[2])]*v[0]*v[1]*u[2];
	if((i1[2]>=0) && (i0[2]<par->nz_here) && (i0[1]>=0) && (i0[1]<par->n_grid) && (i1[0]>=0) && (i1[0]<par->n_grid))
	  mps[id]+=par->grid_vpot[i1[0]+ngx*(i0[1]+par->n_grid*i1[2])]*u[0]*v[1]*u[2];
	if((i1[2]>=0) && (i0[2]<par->nz_here) && (i1[1]>=0) && (i1[1]<par->n_grid) && (i0[0]>=0) && (i0[0]<par->n_grid))
	  mps[id]+=par->grid_vpot[i0[0]+ngx*(i1[1]+par->n_grid*i1[2])]*v[0]*u[1]*u[2];
	if((i1[2]>=0) && (i0[2]<par->nz_here) && (i1[1]>=0) && (i1[1]<par->n_grid) && (i1[0]>=0) && (i1[0]<par->n_grid))
	  mps[id]+=par->grid_vpot[i1[0]+ngx*(i1[1]+par->n_grid*i1[2])]*u[0]*u[1]*u[2];
      }
    } //end omp for
  }
  if(NodeThis==0) timer(2);
}
#endif //_IM_3D22D
#endif //_OLD_IM

void get_IM(ParamCoLoRe *par)
{
  int ipop;

  print_info("*** Getting intensity mapping species\n");
  for(ipop=0;ipop<par->n_im;ipop++)
    get_IM_single(par,ipop);

  printf("Done\n");
  print_info("\n");
}
