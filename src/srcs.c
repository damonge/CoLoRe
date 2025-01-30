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

#ifdef _USE_FAST_LENSING
static int get_r_index_lensing(HealpixShellsAdaptive *sh,double r,int ir_start)
{
  int gotit=0;
  int ir0;
  if(ir_start<0)
    ir0=0;
  else if(ir_start>=sh->nr)
    ir0=sh->nr-1;
  else
    ir0=ir_start;

  while(!gotit) {
    if(ir0==0) {
      if(r<sh->r[1])
        gotit=1;
      else
        ir0++;
    }
    else if(ir0==sh->nr-1) {
      if(r>=sh->r[sh->nr-2]) {
        ir0=sh->nr-2;
        gotit=1;
      }
      else
        ir0--;
    }
    else {
      if(r<sh->r[ir0])
        ir0--;
      else {
        if(r>=sh->r[ir0+1])
          ir0++;
        else
          gotit=1;
      }
    }
  }

  return ir0;
}
#endif //_USE_FAST_LENSING


static inline void cart2sph(double x,double y,double z,double *r,double *cth,double *phi)
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
    *phi=atan2(yn,xn);
    if((*phi)<0)
      (*phi)+=2*M_PI;
  }
}

static double get_rvel(ParamCoLoRe *par,int ix,int iy,int iz,
		       double x0,double y0,double z0,double rr)
{
  double v[3],u[3];
  double idx=par->n_grid/par->l_box;
  long ngx=2*(par->n_grid/2+1);
  long iz_hi=iz+1,iz_lo=iz-1,iz_0=iz;
  long iy_hi=iy+1,iy_lo=iy-1,iy_0=iy;
  long ix_hi=ix+1,ix_lo=ix-1,ix_0=ix;
  if(iy==0) iy_lo=par->n_grid-1;
  if(iy==par->n_grid-1) iy_hi=0;
  if(ix==0) ix_lo=par->n_grid-1;
  if(ix==par->n_grid-1) ix_hi=0;
  iz_0*=ngx*par->n_grid;
  iz_lo*=ngx*par->n_grid;
  iz_hi*=ngx*par->n_grid;
  iy_0*=ngx;
  iy_lo*=ngx;
  iy_hi*=ngx;

  u[0]=x0/rr; u[1]=y0/rr; u[2]=z0/rr;
  v[0]=par->grid_npot[ix_hi+iy_0+iz_0]-par->grid_npot[ix_lo+iy_0+iz_0];
  v[1]=par->grid_npot[ix_0+iy_hi+iz_0]-par->grid_npot[ix_0+iy_lo+iz_0];
  if(iz==0)
    v[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->slice_left[ix_0+iy_0];
  else if(iz==par->nz_here-1)
    v[2]=par->slice_right[ix_0+iy_0]-par->grid_npot[ix_0+iy_0+iz_lo];
  else
    v[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->grid_npot[ix_0+iy_0+iz_lo];

  return 0.5*idx*(v[0]*u[0]+v[1]*u[1]+v[2]*u[2]);
}

static void srcs_set_cartesian_single(ParamCoLoRe *par,int ipop)
{
  int ii,nthr;
  int ngx=2*(par->n_grid/2+1);
  long *np_tot_thr;
  int *nsources=my_calloc(ngx*((long)(par->n_grid*par->nz_here)),sizeof(int));
#ifdef _HAVE_OMP
  nthr=omp_get_max_threads();
#else //_HAVE_OMP
  nthr=1;
#endif //_HAVE_OMP
  np_tot_thr=my_calloc(nthr,sizeof(long));

  print_info(" %d-th galaxy population\n",ipop);
  if(NodeThis==0) timer(0);
  print_info("   Poisson-sampling\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none)			\
  shared(par,np_tot_thr,IThread0,ipop,nthr,nsources)
#endif //_HAVE_OMP
  {
    long iz;
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    double dx=par->l_box/par->n_grid;
    double cell_vol=dx*dx*dx;
    int ngx=2*(par->n_grid/2+1);
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_srcs*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      double z0=(iz+par->iz0_here+0.0)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	double y0=(iy+0.0)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  int npp=0;
	  long index=ix+indexy+indexz;
	  double x0=(ix+0.0)*dx-par->pos_obs[0];
	  double r=sqrt(x0*x0+y0*y0+z0*z0);
	  if(r<par->l_box/2+20.) {
	    double ndens=get_bg(par,r,BG_NZ_SRCS,ipop);
	    if(ndens>0) {
	      double bias=get_bg(par,r,BG_BZ_SRCS,ipop);
	      double dnorm=get_bg(par,r,BG_NORM_SRCS,ipop);
	      double lambda=ndens*cell_vol*bias_model(par->grid_dens[index],bias, get_bg(par, r, BG_TZ_SRCS, ipop))*dnorm;
	      //double lambda=ndens*cell_vol*dnorm;
	      npp=rng_poisson(lambda,rng_thr);
	    }
	  }

	  nsources[index]=npp;
	  np_tot_thr[ithr]+=npp;
	}
      }
    }//end omp for

    end_rng(rng_thr);
  }//end omp parallel

  par->nsources_c_this[ipop]=0;
  for(ii=0;ii<nthr;ii++)
    par->nsources_c_this[ipop]+=np_tot_thr[ii];

  long nsources_total=0;
#ifdef _HAVE_MPI
  MPI_Allreduce(&(par->nsources_c_this[ipop]),&nsources_total,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
#else //_HAVE_MPI
  nsources_total=par->nsources_c_this[ipop];
#endif //_HAVE_MPI

  print_info("   There will be %ld objects in total \n",(long)nsources_total);
#ifdef _DEBUG
  fprintf(par->f_dbg,"MPI task %d has %ld particles\n",NodeThis,(long)(par->nsources_c_this[ipop]));
#endif //_DEBUG

  for(ii=nthr-1;ii>0;ii--) {
    int jj;
    long nh=0;
    for(jj=0;jj<ii;jj++)
      nh+=np_tot_thr[jj];
    np_tot_thr[ii]=nh;
  }
  np_tot_thr[0]=0;
  //np_tot_thr now contains the id of the first particle in the thread

  par->cats_c[ipop]=catalog_cartesian_alloc(par->nsources_c_this[ipop]);

  print_info("   Assigning coordinates\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none)			\
  shared(par,IThread0,np_tot_thr,ipop,nthr,nsources)
#endif //_HAVE_OMP
  {
    long iz;
#ifdef _HAVE_OMP
    int ithr=omp_get_thread_num();
#else //_HAVE_OMP
    int ithr=0;
#endif //_HAVE_OMP
    double dx=par->l_box/par->n_grid;
    int ngx=2*(par->n_grid/2+1);
    unsigned int seed_thr=par->seed_rng+ithr+nthr*(ipop+par->n_srcs*IThread0);
    gsl_rng *rng_thr=init_rng(seed_thr);
    double factor_vel=-par->fgrowth_0/(1.5*par->hubble_0*par->OmegaM);

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      double z0=(iz+par->iz0_here+0.0)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	double y0=(iy+0.0)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  double x0=(ix+0.0)*dx-par->pos_obs[0];
	  long index=ix+indexy+indexz;
	  int npp=nsources[index];
	  if(npp>0) {
	    int ip;
	    double rr=sqrt(x0*x0+y0*y0+z0*z0);
	    double rvel=factor_vel*get_rvel(par,ix,iy,iz,x0,y0,z0,rr);
	    double dz_rsd=rvel*get_bg(par,rr,BG_V1,0);
	    for(ip=0;ip<npp;ip++) {
	      int ax;
	      long pix_id_ring,pix_id_nest;
	      double pos[3];
	      long pid=np_tot_thr[ithr];

	      par->cats_c[ipop]->pos[NPOS_CC*pid+0]=x0+dx*(rng_01(rng_thr)-0.5);
	      par->cats_c[ipop]->pos[NPOS_CC*pid+1]=y0+dx*(rng_01(rng_thr)-0.5);
	      par->cats_c[ipop]->pos[NPOS_CC*pid+2]=z0+dx*(rng_01(rng_thr)-0.5);
	      par->cats_c[ipop]->pos[NPOS_CC*pid+3]=dz_rsd;
	      for(ax=0;ax<3;ax++)
		pos[ax]=par->cats_c[ipop]->pos[NPOS_CC*pid+ax];

	      vec2pix_ring(par->nside_base,pos,&pix_id_ring);
	      ring2nest(par->nside_base,pix_id_ring,&pix_id_nest);
	      par->cats_c[ipop]->ipix[pid]=pix_id_nest;
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
  free(nsources);
}

void srcs_set_cartesian(ParamCoLoRe *par)
{
  int ipop;

  //First, compute lensing Hessian
  print_info("*** Getting point sources\n");
  for(ipop=0;ipop<par->n_srcs;ipop++)
    srcs_set_cartesian_single(par,ipop);
  print_info("\n");
}

static void srcs_distribute_single(ParamCoLoRe *par,int ipop)
{
  int ii;
  long *ns_to_nodes=my_calloc(NNodes,sizeof(long));
  long *i0_in_nodes=my_calloc(NNodes,sizeof(long));
  long *ns_transfer_matrix=my_calloc(NNodes*NNodes,sizeof(long));
  CatalogCartesian *cat=par->cats_c[ipop];

  //Count how many objects need to be sent to each node
  for(ii=0;ii<par->nsources_c_this[ipop];ii++) {
    int i_node=cat->ipix[ii]%NNodes;
    ns_to_nodes[i_node]++;
  }

  //Gather all transfers into a 2D array
#ifdef _HAVE_MPI
  MPI_Allgather(ns_to_nodes,NNodes,MPI_LONG,
		ns_transfer_matrix,NNodes,MPI_LONG,MPI_COMM_WORLD);
#else //_HAVE_MPI
  ns_transfer_matrix[0]=ns_to_nodes[0];
#endif //_HAVE_MPI

  //Calculate position of each particle going to each node
  for(ii=0;ii<NNodes;ii++) {
    int jj;
    i0_in_nodes[ii]=0;
    for(jj=0;jj<ii;jj++)
      i0_in_nodes[ii]+=ns_to_nodes[jj];
  }
  for(ii=0;ii<NNodes;ii++)
    ns_to_nodes[ii]=i0_in_nodes[ii];

  //Reorganize positions so they are ordered by destination node
  float *pos_ordered=my_malloc(NPOS_CC*par->nsources_c_this[ipop]*sizeof(float));
  for(ii=0;ii<par->nsources_c_this[ipop];ii++) {
    int ax;
    int i_node=cat->ipix[ii]%NNodes;
    long id_here=ns_to_nodes[i_node];
    for(ax=0;ax<NPOS_CC;ax++)
      pos_ordered[NPOS_CC*id_here+ax]=cat->pos[NPOS_CC*ii+ax];
    ns_to_nodes[i_node]++;
  }

  //Count particles this node will receive
  par->nsources_this[ipop]=0;
  for(ii=0;ii<NNodes;ii++)
    par->nsources_this[ipop]+=ns_transfer_matrix[ii*NNodes+NodeThis];

  //Loop through all nodes and receive particles
  //Free up old catalog and create new one
  catalog_cartesian_free(par->cats_c[ipop]);
  par->cats_c[ipop]=catalog_cartesian_alloc(par->nsources_this[ipop]);
  par->nsources_c_this[ipop]=par->nsources_this[ipop];
  long i_sofar=0;
  for(ii=0;ii<NNodes;ii++) {
    int node_to  =(NodeRight+ii+NNodes)%NNodes;
    int node_from=(NodeLeft -ii+NNodes)%NNodes;
    int npart_send=ns_transfer_matrix[NodeThis*NNodes+node_to];
    int npart_recv=ns_transfer_matrix[node_from*NNodes+NodeThis];
    //    print_info("MPI task %d: %d-th iteration. to->%d from->%d.",NodeThis,ii,node_to,node_from);
    //    print_info(" Should get %07ld objects, and will send %07ld.\n",npart_recv,npart_send);
#ifdef _HAVE_MPI
    MPI_Sendrecv(&(pos_ordered[NPOS_CC*i0_in_nodes[node_to]]),NPOS_CC*npart_send,MPI_FLOAT,node_to  ,ii,
		 &(par->cats_c[ipop]->pos[NPOS_CC*i_sofar])  ,NPOS_CC*npart_recv,MPI_FLOAT,node_from,ii,
		 MPI_COMM_WORLD,MPI_STATUS_IGNORE);
#else //_HAVE_MPI
    memcpy(&(par->cats_c[ipop]->pos[NPOS_CC*i_sofar]),&(pos_ordered[NPOS_CC*i0_in_nodes[node_to]]),
	   NPOS_CC*npart_send*sizeof(float));
#endif //_HAVE_MPI
    i_sofar+=npart_recv;
    //    print_info(" %07ld sources gathered so far\n",i_sofar);
  }

  free(pos_ordered);
  free(i0_in_nodes);
  free(ns_to_nodes);
  free(ns_transfer_matrix);
}

void srcs_distribute(ParamCoLoRe *par)
{
  if(NodeThis==0) timer(0);
  int ipop;
  print_info("*** Re-distributing sources across nodes\n");
  for(ipop=0;ipop<par->n_srcs;ipop++)
    srcs_distribute_single(par,ipop);
  if(NodeThis==0) timer(2);
  print_info("\n");
}

static void srcs_get_local_properties_single(ParamCoLoRe *par,int ipop)
{
  par->cats[ipop]=catalog_alloc(par->cats_c[ipop]->nsrc,
                                par->lensing_srcs[ipop],
                                par->skw_srcs[ipop],
                                par->skw_gauss[ipop],
				par->r_max,par->n_grid);

#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par,ipop)
#endif //_HAVE_OMP
  {
    int ii;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<par->cats[ipop]->nsrc;ii++) {
      double r,cth,phi;
      float *pos=&(par->cats_c[ipop]->pos[NPOS_CC*ii]);
      cart2sph(pos[0],pos[1],pos[2],&r,&cth,&phi);
      par->cats[ipop]->srcs[ii].ra=RTOD*phi;
      par->cats[ipop]->srcs[ii].dec=90-RTOD*acos(cth);
      par->cats[ipop]->srcs[ii].z0=get_bg(par,r,BG_Z,0);
      par->cats[ipop]->srcs[ii].dz_rsd=pos[3];
      par->cats[ipop]->srcs[ii].e1=-1;
      par->cats[ipop]->srcs[ii].e2=-1;
    }//end omp for
  }//end omp parallel
}

void srcs_get_local_properties(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_srcs;ipop++)
    srcs_get_local_properties_single(par,ipop);
}

static void srcs_beams_preproc_single(ParamCoLoRe *par,int ipop)
{
#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par,ipop)
#endif //_HAVE_OMP
  {
    int ii;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<par->cats[ipop]->nsrc;ii++) {
      par->cats[ipop]->srcs[ii].dz_rsd=0;
      par->cats[ipop]->srcs[ii].e1=0;
      par->cats[ipop]->srcs[ii].e2=0;
    }//end omp for
  }//end omp parallel
}

void srcs_beams_preproc(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_srcs;ipop++)
    srcs_beams_preproc_single(par,ipop);
}

static void srcs_get_beam_properties_single(ParamCoLoRe *par,int ipop)
{
  Catalog *cat=par->cats[ipop];
  CatalogCartesian *catc=par->cats_c[ipop];

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,cat,catc)
#endif //_HAVE_OMP
  {
    long ip;
    double idx=par->n_grid/par->l_box;

    //Kernels for the LOS integrals
#ifndef _USE_FAST_LENSING
    double *fac_r_0=NULL,*fac_r_1=NULL,*fac_r_2=NULL;
    if(cat->has_lensing) {
      fac_r_0=my_malloc(cat->nr*sizeof(double));
      fac_r_1=my_malloc(cat->nr*sizeof(double));
      fac_r_2=my_malloc(cat->nr*sizeof(double));

      for(ip=0;ip<cat->nr;ip++) {
	double rm=(ip+0.5)*cat->dr;
	double pg=get_bg(par,rm,BG_D1,0)*(1+get_bg(par,rm,BG_Z,0));
	fac_r_0[ip]=pg*cat->dr;
	fac_r_1[ip]=rm*pg*cat->dr;
	fac_r_2[ip]=rm*rm*pg*cat->dr;
      }
    }
#endif //_USE_FAST_LENSING

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ip=0;ip<catc->nsrc;ip++) {
      int ax,added;
      double xn[3],u[3];
      flouble v[3],vr,dens,gauss;
      float *pos=&(catc->pos[NPOS_CC*ip]);
      double r=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
      double ir=1./(MAX(r,0.001));

      for(ax=0;ax<3;ax++) {//TODO these could be precomputed in the cartesian catalog
	xn[ax]=(pos[ax]+par->pos_obs[ax])*idx;
	u[ax]=pos[ax]*ir;
      }

      //Compute RSD
      added=interpolate_from_grid(par,xn,NULL,v,NULL,NULL,NULL,RETURN_VEL,INTERP_TYPE_SKW);
      if(added) {
	vr=0.5*idx*(v[0]*u[0]+v[1]*u[1]+v[2]*u[2]);
	cat->srcs[ip].dz_rsd+=vr;
      }

      //Fill up skewers
      if(cat->has_skw) {
      	int i_r,i_r_max=MIN((int)(r*cat->idr+0.5),cat->nr-1);
      	long offp=ip*cat->nr;
      	for(i_r=0;i_r<=i_r_max;i_r++) {
      	  double rm=(i_r+0.5)*cat->dr;
          for(ax=0;ax<3;ax++)
      	    xn[ax]=(rm*u[ax]+par->pos_obs[ax])*idx;
	  if(cat->skw_gauss) {
	    added=interpolate_from_grid(par,xn,&dens,v,NULL,NULL,&gauss,RETURN_GAUSS | RETURN_VEL,INTERP_TYPE_SKW);
	  }
	  else {
	    added=interpolate_from_grid(par,xn,&dens,v,NULL,NULL,&gauss,RETURN_DENS | RETURN_VEL,INTERP_TYPE_SKW);
	  }
      	  if(added) {
      	    vr=0.5*idx*(v[0]*u[0]+v[1]*u[1]+v[2]*u[2]);
            if(cat->skw_gauss)
              cat->g_skw[offp+i_r]+=gauss;
            else
              cat->d_skw[offp+i_r]+=dens;
      	    cat->v_skw[offp+i_r]+=vr;
      	  }
      	}
      }

#ifndef _USE_FAST_LENSING
      //Compute lensing
      if(cat->has_lensing) {
	//Compute linear transformations needed for lensing
        double u_x[3], u_y[3];
	double r_k[6], r_e1[6],r_e2[6];
	double cth_h=1,sth_h=0,cph_h=1,sph_h=0;
	double prefac=idx*idx*ir;  // 1 / (Dx^2 * r)
        double prefac_m=0.5*idx*ir;  // 1 / (2 * Dx * r)

	cth_h=u[2];
	if(cth_h>=1) cth_h=1;
	if(cth_h<=-1) cth_h=-1;
	sth_h=sqrt((1-cth_h)*(1+cth_h));
	if(sth_h!=0) {
	  cph_h=u[0]/sth_h;
	  sph_h=u[1]/sth_h;
	}

        u_x[0]=cth_h*cph_h*prefac_m;
        u_x[1]=cth_h*sph_h*prefac_m;
        u_x[2]=-sth_h*prefac_m;

        u_y[0]=-sph_h*prefac_m;
        u_y[1]=cph_h*prefac_m;
        u_y[2]=0;

        r_k[0] =(cth_h*cth_h*cph_h*cph_h+sph_h*sph_h)*prefac;
        r_k[1] =(2*cph_h*sph_h*(cth_h*cth_h-1))*prefac;
        r_k[2] =(-2*cth_h*sth_h*cph_h)*prefac;
        r_k[3] =(cth_h*cth_h*sph_h*sph_h+cph_h*cph_h)*prefac;
        r_k[4] =(-2*cth_h*sth_h*sph_h)*prefac;
        r_k[5] =(sth_h*sth_h)*prefac;

	r_e1[0]=(cth_h*cth_h*cph_h*cph_h-sph_h*sph_h)*prefac;
	r_e1[1]=(2*cph_h*sph_h*(cth_h*cth_h+1))*prefac;
	r_e1[2]=(-2*cth_h*sth_h*cph_h)*prefac;
	r_e1[3]=(cth_h*cth_h*sph_h*sph_h-cph_h*cph_h)*prefac;
	r_e1[4]=(-2*cth_h*sth_h*sph_h)*prefac;
	r_e1[5]=(sth_h*sth_h)*prefac;

	r_e2[0]=(-2*cth_h*cph_h*sph_h)*prefac;
	r_e2[1]=(2*cth_h*(cph_h*cph_h-sph_h*sph_h))*prefac;
	r_e2[2]=(2*sth_h*sph_h)*prefac;
	r_e2[3]=(2*cth_h*sph_h*cph_h)*prefac;
	r_e2[4]=(-2*sth_h*cph_h)*prefac;
	r_e2[5]=0;

	//Integrate along the LOS
	flouble tp[6],vp[3];
	double dtx=0,dty=0,kp=0,e1=0,e2=0;
	int i_r,i_r_max=MIN((int)(r*cat->idr+0.5),cat->nr-1);
	for(i_r=0;i_r<=i_r_max;i_r++) {
	  double rm=(i_r+0.5)*cat->dr;
	  for(ax=0;ax<3;ax++)
	    xn[ax]=(rm*u[ax]+par->pos_obs[ax])*idx;
	  added=interpolate_from_grid(par,xn,NULL,vp,tp,NULL,NULL,RETURN_VEL | RETURN_TID,
                                      INTERP_TYPE_LENSING);
	  if(added) {
	    double fr=fac_r_1[i_r]*r-fac_r_2[i_r];
	    double frm=2*(fac_r_0[i_r]*r-fac_r_1[i_r]);
	    double dotvx=0,dotvy=0,dotk=0,dote1=0,dote2=0;
	    for(ax=0;ax<6;ax++) {
	      dote1+=r_e1[ax]*tp[ax];
	      dote2+=r_e2[ax]*tp[ax];
	      dotk+=r_k[ax]*tp[ax];
	    }
	    for(ax=0;ax<3;ax++) {
              dotvx+=u_x[ax]*vp[ax];
              dotvy+=u_y[ax]*vp[ax];
            }
	    e1+=dote1*fr;
	    e2+=dote2*fr;
	    kp+=dotk*fr;
            dtx+=dotvx*frm;
            dty+=dotvy*frm;
	  }
	}
	cat->srcs[ip].e1+=e1;
	cat->srcs[ip].e2+=e2;
	cat->srcs[ip].kappa+=kp;
        cat->srcs[ip].dra+=dty;
        cat->srcs[ip].ddec+=dtx;
      }
#endif //_USE_FAST_LENSING
    }//end omp for
#ifndef _USE_FAST_LENSING
    if(cat->has_lensing) {
      free(fac_r_0);
      free(fac_r_1);
      free(fac_r_2);
    }
#endif //_USE_FAST_LENSING
  }//end omp parallel
}

void srcs_get_beam_properties(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_srcs;ipop++)
    srcs_get_beam_properties_single(par,ipop);
}

static void srcs_beams_postproc_single(ParamCoLoRe *par,int ipop)
{
  Catalog *cat=par->cats[ipop];
  CatalogCartesian *catc=par->cats_c[ipop];

#ifdef _HAVE_OMP
#pragma omp parallel default(none)              \
  shared(par,cat,catc,NNodes,NodeThis)
#endif //_HAVE_OMP
  {
    int ii;
    double factor_vel=-par->fgrowth_0/(1.5*par->hubble_0*par->OmegaM);
#ifdef _USE_FAST_LENSING
    int ir_s=0;
    HealpixShellsAdaptive *smap=NULL;
    if(cat->has_lensing)
      smap = par->smap;
#endif //_USE_FAST_LENSING

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<cat->nsrc;ii++) {
      double z=cat->srcs[ii].z0;
      double r=r_of_z(par,z);
      double vg=get_bg(par,r,BG_V1,0);

      //RSDs
      cat->srcs[ii].dz_rsd*=vg*factor_vel;

      //Lensing
      if(cat->has_lensing) {
#ifdef _USE_FAST_LENSING
        int ax;
        long ipix_up,ipix_lo,ibase,ibase_here;
        double u[3];
        double g1_lo,g1_up,g2_lo,g2_up;
        double kp_lo,kp_up;
        double dx_lo,dx_up,dy_lo,dy_up;
        //Find lower shell index
        ir_s=get_r_index_lensing(smap, r, ir_s);
        if(ir_s>=smap->nr-1)
          report_error(1, "SHIT\n");

        //Find intervals
        double h = (r - smap->r[ir_s]) / (smap->r[ir_s+1] - smap->r[ir_s]);

        //Find pixel index
        //Find base this galaxy belongs to
        for(ax=0;ax<3;ax++)
          u[ax]=catc->pos[NPOS_CC*ii+ax];
        vec2pix_nest(par->nside_base,u,&ibase);
        if(ibase%NNodes!=NodeThis)
          report_error(1,"Bad base!!\n");
        //Find base index within this node
        ibase_here=(ibase-NodeThis)/NNodes;

        //Find pixel in lower edge
        vec2pix_nest(smap->nside[ir_s],u,&ipix_lo);
        //Offset to pixel indices in this node
        ipix_lo-=ibase*smap->num_pix_per_beam[ir_s];
	if((ipix_lo<0) || (ipix_lo>=smap->num_pix_per_beam[ir_s]))
	  report_error(1,"Bad pixel!!\n");
        //Find lensing at edges
        g1_lo=smap->data[ibase_here][ir_s][5*ipix_lo+0];
        g2_lo=smap->data[ibase_here][ir_s][5*ipix_lo+1];
        kp_lo=smap->data[ibase_here][ir_s][5*ipix_lo+2];
        dx_lo=smap->data[ibase_here][ir_s][5*ipix_lo+3];
        dy_lo=smap->data[ibase_here][ir_s][5*ipix_lo+4];

        //Same in upper edge
        vec2pix_nest(smap->nside[ir_s+1],u,&ipix_up);
        ipix_up-=ibase*smap->num_pix_per_beam[ir_s+1];
	if((ipix_up<0) || (ipix_up>=smap->num_pix_per_beam[ir_s+1]))
	  report_error(1,"Bad pixel!!\n");
        g1_up=smap->data[ibase_here][ir_s+1][2*ipix_up+0];
        g2_up=smap->data[ibase_here][ir_s+1][2*ipix_up+1];
        kp_up=smap->data[ibase_here][ir_s+1][2*ipix_up+2];
        dx_up=smap->data[ibase_here][ir_s+1][2*ipix_up+3];
        dy_up=smap->data[ibase_here][ir_s+1][2*ipix_up+4];

        //Interpolate
        cat->srcs[ii].e1=g1_lo*(1-h)+g1_up*h;
        cat->srcs[ii].e2=g2_lo*(1-h)+g2_up*h;
        cat->srcs[ii].kappa=kp_lo*(1-h)+kp_up*h;
        cat->srcs[ii].dra=dy_lo*(1-h)+dy_up*h;
        cat->srcs[ii].ddec=dx_lo*(1-h)+dx_up*h;
#endif //_USE_FAST_LENSING
        cat->srcs[ii].dra*=RTOD;
        cat->srcs[ii].ddec*=RTOD;
      }

      //Skewers
      if(cat->has_skw) {
	int i_r,i_r_max=MAX((int)(r*cat->idr+0.5),cat->nr-1);
	long offp=ii*cat->nr;
	for(i_r=0;i_r<=i_r_max;i_r++) {
	  vg=get_bg(par,(i_r+0.5)*cat->dr,BG_V1,0);
	  cat->v_skw[offp+i_r]*=vg*factor_vel;
	}
      }
    }//end omp for
  }//end omp parallel
}

void srcs_beams_postproc(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_srcs;ipop++)
    srcs_beams_postproc_single(par,ipop);
}
