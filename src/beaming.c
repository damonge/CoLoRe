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

#define IND_XX 0
#define IND_XY 1
#define IND_XZ 2
#define IND_YY 3
#define IND_YZ 4
#define IND_ZZ 5

static void get_element(ParamCoLoRe *par,long ix,long iy,long iz,
      flouble *d,flouble v[3],flouble dfl[3],flouble t[6],flouble *pdot,flouble *g,
      int flag_return)
{
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

  //Get density
  if(flag_return & RETURN_DENS)
    *d=par->grid_dens[ix_0+iy_0+iz_0];

  //Get gaussian
  if(flag_return & RETURN_GAUSS) {
    if (par->dens_type!=DENS_TYPE_LGNR)
      report_error(1,"Cannot write Gaussian skewers with density type %d\n",par->dens_type);	    
    //Get D and r
    flouble dx=par->l_box/par->n_grid;
    flouble z0=(iz+par->iz0_here+0.0)*dx-par->pos_obs[2];
    flouble y0=(iy+0.0)*dx-par->pos_obs[1];
    flouble x0=(ix+0.0)*dx-par->pos_obs[0];
    double r=sqrt(x0*x0+y0*y0+z0*z0);
    double dg=get_bg(par,r,BG_D1,0);
    *g=(log(1+par->grid_dens[ix_0+iy_0+iz_0]))/dg+dg*(par->sigma2_gauss)*0.5;
  }

  //Get velocity
  if(flag_return & RETURN_VEL) {
    if (par->lpt_vels){
       v[0]=par->grid_velx[ix_0+iy_0+iz_0];
       v[1]=par->grid_vely[ix_0+iy_0+iz_0];
       v[2]=par->grid_velz[ix_0+iy_0+iz_0];
    } else {
       v[0]=par->grid_npot[ix_hi+iy_0+iz_0]-par->grid_npot[ix_lo+iy_0+iz_0];
       v[1]=par->grid_npot[ix_0+iy_hi+iz_0]-par->grid_npot[ix_0+iy_lo+iz_0];
       if(iz==0)
          v[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->slice_left[ix_0+iy_0];
       else if(iz==par->nz_here-1)
         v[2]=par->slice_right[ix_0+iy_0]-par->grid_npot[ix_0+iy_0+iz_lo];
       else
         v[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->grid_npot[ix_0+iy_0+iz_lo];
    }
  }

  //Get deflection
  if(flag_return & RETURN_DFL) {
    dfl[0]=par->grid_npot[ix_hi+iy_0+iz_0]-par->grid_npot[ix_lo+iy_0+iz_0];
    dfl[1]=par->grid_npot[ix_0+iy_hi+iz_0]-par->grid_npot[ix_0+iy_lo+iz_0];
    if(iz==0)
      dfl[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->slice_left[ix_0+iy_0];
    else if(iz==par->nz_here-1)
      dfl[2]=par->slice_right[ix_0+iy_0]-par->grid_npot[ix_0+iy_0+iz_lo];
    else
      dfl[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->grid_npot[ix_0+iy_0+iz_lo];
  }

  //Get ISW
  if(flag_return & RETURN_PDOT)
    *pdot=par->grid_npot[ix_0+iy_0+iz_0];

  //Get tidal tensor
  if(flag_return & RETURN_TID) {
    t[IND_XX]=(par->grid_npot[ix_hi+iy_0+iz_0]+par->grid_npot[ix_lo+iy_0+iz_0]-
              2*par->grid_npot[ix_0+iy_0+iz_0]);
    t[IND_YY]=(par->grid_npot[ix_0+iy_hi+iz_0]+par->grid_npot[ix_0+iy_lo+iz_0]-
              2*par->grid_npot[ix_0+iy_0+iz_0]);
    t[IND_XY]=0.25*(par->grid_npot[ix_hi+iy_hi+iz_0]+par->grid_npot[ix_lo+iy_lo+iz_0]-
                   par->grid_npot[ix_hi+iy_lo+iz_0]-par->grid_npot[ix_lo+iy_hi+iz_0]);
    if(iz==0) {
      t[IND_ZZ]=(par->grid_npot[ix_0+iy_0+iz_hi]+par->slice_left[ix_0+iy_0]-
                2*par->grid_npot[ix_0+iy_0+iz_0]);
      t[IND_XZ]=0.25*(par->grid_npot[ix_hi+iy_0+iz_hi]+par->slice_left[ix_lo+iy_0]-
                     par->slice_left[ix_hi+iy_0]-par->grid_npot[ix_lo+iy_0+iz_hi]);
      t[IND_YZ]=0.25*(par->grid_npot[ix_0+iy_hi+iz_hi]+par->slice_left[ix_0+iy_lo]-
                     par->slice_left[ix_0+iy_hi]-par->grid_npot[ix_0+iy_lo+iz_hi]);
    }
    else if(iz==par->nz_here-1) {
      t[IND_ZZ]=(par->slice_right[ix_0+iy_0]+par->grid_npot[ix_0+iy_0+iz_lo]-
                2*par->grid_npot[ix_0+iy_0+iz_0]);
      t[IND_XZ]=0.25*(par->slice_right[ix_hi+iy_0]+par->grid_npot[ix_lo+iy_0+iz_lo]-
                     par->grid_npot[ix_hi+iy_0+iz_lo]-par->slice_right[ix_lo+iy_0]);
      t[IND_YZ]=0.25*(par->slice_right[ix_0+iy_hi]+par->grid_npot[ix_0+iy_lo+iz_lo]-
                     par->grid_npot[ix_0+iy_hi+iz_lo]-par->slice_right[ix_0+iy_lo]);
    }
    else {
      t[IND_ZZ]=(par->grid_npot[ix_0+iy_0+iz_hi]+par->grid_npot[ix_0+iy_0+iz_lo]-
                2*par->grid_npot[ix_0+iy_0+iz_0]);
      t[IND_XZ]=0.25*(par->grid_npot[ix_hi+iy_0+iz_hi]+par->grid_npot[ix_lo+iy_0+iz_lo]-
                     par->grid_npot[ix_hi+iy_0+iz_lo]-par->grid_npot[ix_lo+iy_0+iz_hi]);
      t[IND_YZ]=0.25*(par->grid_npot[ix_0+iy_hi+iz_hi]+par->grid_npot[ix_0+iy_lo+iz_lo]-
                     par->grid_npot[ix_0+iy_hi+iz_lo]-par->grid_npot[ix_0+iy_lo+iz_hi]);
    }
  }
}

//WARNING!!!! !: x should go from 0 to ngrid-1!!!!! Remove this when taken care of
int interpolate_from_grid(ParamCoLoRe *par,double *x,
                          flouble *d,flouble v[3],flouble dfl[3],
			  flouble t[6],flouble *pd,flouble *g,
                          int flag_return,int interp_type)
{
  long ix0[3];
  double h0x[3];
  int ax,added_anything=0;

  //Initialize output
  if(flag_return & RETURN_DENS)
    *d=0;

  if(flag_return & RETURN_GAUSS)
    *g=0;

  if(flag_return & RETURN_VEL) {
    for(ax=0;ax<3;ax++)
      v[ax]=0;
  }

  if(flag_return & RETURN_DFL) {
    for(ax=0;ax<3;ax++)
      dfl[ax]=0;
  }

  if(flag_return & RETURN_TID) {
    for(ax=0;ax<6;ax++)
      t[ax]=0;
  }

  if(flag_return & RETURN_PDOT)
    *pd=0;

  if(interp_type==INTERP_NGP) {
    for(ax=0;ax<3;ax++) {
      ix0[ax]=(long)(x[ax]+0.5);
      if(ix0[ax]>=par->n_grid)
        ix0[ax]-=par->n_grid;
      else if(ix0[ax]<0)
        ix0[ax]+=par->n_grid;
      h0x[ax]=1.;
    }
    ix0[2]-=par->iz0_here;

    if((ix0[2]>=0) && (ix0[2]<par->nz_here)) {
      flouble d_000,g_000,v_000[3],dfl_000[3],t_000[6],pd_000;
      flouble w_000=h0x[2]*h0x[1]*h0x[0];

      added_anything=1;
      get_element(par,ix0[0],ix0[1],ix0[2],&d_000,v_000,dfl_000,t_000,&pd_000,&g_000,flag_return);
      if(flag_return & RETURN_DENS)
        *d+=d_000*w_000;

      if(flag_return & RETURN_GAUSS)
        *g+=g_000*w_000;

      if(flag_return & RETURN_VEL) {
        for(ax=0;ax<3;ax++)
          v[ax]+=v_000[ax]*w_000;
      }
      if(flag_return & RETURN_DFL) {
        for(ax=0;ax<3;ax++)
          dfl[ax]+=dfl_000[ax]*w_000;
      }
      if(flag_return & RETURN_TID) {
        for(ax=0;ax<6;ax++)
          t[ax]+=t_000[ax]*w_000;
      }
      if(flag_return & RETURN_PDOT)
        *pd+=pd_000*w_000;
    }
  }
  else {
    long ix1[3];
    flouble h1x[3];

    //Trilinear interpolation
    for(ax=0;ax<3;ax++) {
      ix0[ax]=(long)(x[ax]);
      h0x[ax]=x[ax]-ix0[ax];
      h1x[ax]=1-h0x[ax];
      ix1[ax]=ix0[ax]+1;
      if(ix0[ax]>=par->n_grid)
        ix0[ax]-=par->n_grid;
      else if(ix0[ax]<0)
        ix0[ax]+=par->n_grid;
      if(ix1[ax]>=par->n_grid)
        ix1[ax]-=par->n_grid;
      else if(ix1[ax]<0)
        ix1[ax]+=par->n_grid;
    }
    ix0[2]-=par->iz0_here;
    ix1[2]-=par->iz0_here;

    if((ix0[2]>=0) && (ix0[2]<par->nz_here)) {
      flouble d_000,v_000[3],dfl_000[3],t_000[6],pd_000,g_000;
      flouble d_001,v_001[3],dfl_001[3],t_001[6],pd_001,g_001;
      flouble d_010,v_010[3],dfl_010[3],t_010[6],pd_010,g_010;
      flouble d_011,v_011[3],dfl_011[3],t_011[6],pd_011,g_011;
      flouble w_000=h1x[2]*h1x[1]*h1x[0];
      flouble w_001=h1x[2]*h1x[1]*h0x[0];
      flouble w_010=h1x[2]*h0x[1]*h1x[0];
      flouble w_011=h1x[2]*h0x[1]*h0x[0];

      added_anything=1;
      get_element(par,ix0[0],ix0[1],ix0[2],&d_000,v_000,dfl_000,t_000,&pd_000,&g_000,flag_return);
      get_element(par,ix1[0],ix0[1],ix0[2],&d_001,v_001,dfl_001,t_001,&pd_001,&g_001,flag_return);
      get_element(par,ix0[0],ix1[1],ix0[2],&d_010,v_010,dfl_010,t_010,&pd_010,&g_010,flag_return);
      get_element(par,ix1[0],ix1[1],ix0[2],&d_011,v_011,dfl_011,t_011,&pd_011,&g_011,flag_return);
      if(flag_return & RETURN_DENS)
        *d+=(d_000*w_000+d_001*w_001+d_010*w_010+d_011*w_011);
      if(flag_return & RETURN_GAUSS)
        *g+=(g_000*w_000+g_001*w_001+g_010*w_010+g_011*w_011);
      if(flag_return & RETURN_VEL) {
        for(ax=0;ax<3;ax++)
          v[ax]+=(v_000[ax]*w_000+v_001[ax]*w_001+v_010[ax]*w_010+v_011[ax]*w_011);
      }
      if(flag_return & RETURN_DFL) {
        for(ax=0;ax<3;ax++)
          dfl[ax]+=(dfl_000[ax]*w_000+dfl_001[ax]*w_001+dfl_010[ax]*w_010+dfl_011[ax]*w_011);
      }
      if(flag_return & RETURN_TID) {
        for(ax=0;ax<6;ax++)
          t[ax]+=(t_000[ax]*w_000+t_001[ax]*w_001+t_010[ax]*w_010+t_011[ax]*w_011);
      }
      if(flag_return & RETURN_PDOT)
        *pd+=(pd_000*w_000+pd_001*w_001+pd_010*w_010+pd_011*w_011);
    }
    if((ix1[2]>=0) && (ix1[2]<par->nz_here)) {
      flouble d_100,v_100[3],dfl_100[3],t_100[6],pd_100,g_100;
      flouble d_101,v_101[3],dfl_101[3],t_101[6],pd_101,g_101;
      flouble d_110,v_110[3],dfl_110[3],t_110[6],pd_110,g_110;
      flouble d_111,v_111[3],dfl_111[3],t_111[6],pd_111,g_111;
      flouble w_100=h0x[2]*h1x[1]*h1x[0];
      flouble w_101=h0x[2]*h1x[1]*h0x[0];
      flouble w_110=h0x[2]*h0x[1]*h1x[0];
      flouble w_111=h0x[2]*h0x[1]*h0x[0];

      added_anything=1;
      get_element(par,ix0[0],ix0[1],ix1[2],&d_100,v_100,dfl_100,t_100,&pd_100,&g_100,flag_return);
      get_element(par,ix1[0],ix0[1],ix1[2],&d_101,v_101,dfl_101,t_101,&pd_101,&g_101,flag_return);
      get_element(par,ix0[0],ix1[1],ix1[2],&d_110,v_110,dfl_110,t_110,&pd_110,&g_110,flag_return);
      get_element(par,ix1[0],ix1[1],ix1[2],&d_111,v_111,dfl_111,t_111,&pd_111,&g_111,flag_return);
      if(flag_return & RETURN_DENS)
        *d+=(d_100*w_100+d_101*w_101+d_110*w_110+d_111*w_111);
      if(flag_return & RETURN_GAUSS)
        *g+=(g_100*w_100+g_101*w_101+g_110*w_110+g_111*w_111);
      if(flag_return & RETURN_VEL) {
        for(ax=0;ax<3;ax++)
          v[ax]+=(v_100[ax]*w_100+v_101[ax]*w_101+v_110[ax]*w_110+v_111[ax]*w_111);
      }
      if(flag_return & RETURN_DFL) {
        for(ax=0;ax<3;ax++)
          dfl[ax]+=(dfl_100[ax]*w_100+dfl_101[ax]*w_101+dfl_110[ax]*w_110+dfl_111[ax]*w_111);
      }
      if(flag_return & RETURN_TID) {
        for(ax=0;ax<6;ax++)
          t[ax]+=(t_100[ax]*w_100+t_101[ax]*w_101+t_110[ax]*w_110+t_111[ax]*w_111);
      }
      if(flag_return & RETURN_PDOT)
        *pd+=(pd_100*w_100+pd_101*w_101+pd_110*w_110+pd_111*w_111);
    }
  }

  return added_anything;
}

//All the beam operations on the slab that the task holds now (par->nz_here, par->iz0_here and the fields)
static void beam_operations_one_slab(ParamCoLoRe *par)
{
  if(par->do_kappa)
    kappa_get_beam_properties(par);
#ifdef _USE_FAST_LENSING
  if(par->do_lensing)
    lensing_get_beam_properties(par);
#endif //_USE_FAST_LENSING
  if(par->do_isw)
    isw_get_beam_properties(par);
  if(par->do_srcs)
    srcs_get_beam_properties(par);
  if(par->do_imap)
    imap_get_beam_properties(par);
  if(par->do_cstm)
    cstm_get_beam_properties(par);
}

#ifdef _HAVE_MPI
static double t_ring_mpi=0,t_ring_copy=0; //Time spent in the ring exchange (for COLORE_TIMING)

//Copy of a slab after the exchange, shared among the OpenMP threads (it was a single memcpy, which
//runs at the speed of one core and is done 5 times per round on slabs of GBs)
static void copy_slab(flouble *dest,flouble *src,long count)
{
#define COPY_CHUNK 1048576
  long nchunks=(count+COPY_CHUNK-1)/COPY_CHUNK,ic;
#ifdef _HAVE_OMP
#pragma omp parallel for schedule(static)
#endif //_HAVE_OMP
  for(ic=0;ic<nchunks;ic++) {
    long lo=ic*COPY_CHUNK,hi=lo+COPY_CHUNK;
    if(hi>count) hi=count;
    memcpy(&(dest[lo]),&(src[lo]),(hi-lo)*sizeof(flouble));
  }
}

static void mpi_sendrecv_wrap(flouble *data,flouble *buff,long count,int tag)
{
  // still need to compile even if never called
#define SENDRECV_BATCH 1073741824
  int remainder;
  long i_sofar=0;
  double t0=MPI_Wtime();
  while(i_sofar+SENDRECV_BATCH<count) {
    MPI_Sendrecv(&(data[i_sofar]),SENDRECV_BATCH,FLOUBLE_MPI,NodeRight,tag,
                &(buff[i_sofar]),SENDRECV_BATCH,FLOUBLE_MPI,NodeLeft ,tag,
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    i_sofar+=SENDRECV_BATCH;
  }
  remainder=(int)(count-i_sofar);
  if(remainder>0) {
    MPI_Sendrecv(&(data[i_sofar]),remainder,FLOUBLE_MPI,NodeRight,tag,
                &(buff[i_sofar]),remainder,FLOUBLE_MPI,NodeLeft ,tag,
                MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  double t1=MPI_Wtime();
  copy_slab(data,buff,count);
  t_ring_mpi+=t1-t0;
  t_ring_copy+=MPI_Wtime()-t1;
}
#endif //_HAVE_MPI

#ifdef _HAVE_MPI
//Ring exchange of the slabs in which the transfer of the next slab runs while the current one is
//processed, and no copy back: two sets of buffers, swapped after every round.
//Only the fields that some active beam operation reads go around the ring (see ring_fields_needed):
//the density, the velocities if they come from the LPT, and the Newtonian potential (its buffer has the
//two neighbouring planes at the end, which are used by the tidal tensor, the velocities and the ISW
//term, so they travel with the slab). The pointers of the fields that are not sent are NULL while the
//slabs are processed, so that a read of one of them stops the run instead of using the wrong slab.
//The own slab is processed first, then the slab of the left neighbour, and so on. Every beam operation
//adds the contribution of each slab to what it has accumulated, so the order of the slabs only changes
//the rounding of the sums (for the source skewers, with at most two slabs per point, not even that).
static void ring_fields_needed(ParamCoLoRe *par,int need_npot,int *send_dens,int *send_vel,int *send_npot)
{
  int ipop,srcs_lensing=0;
  if(par->do_srcs) {
    for(ipop=0;ipop<par->n_srcs;ipop++) {
      if(par->lensing_srcs[ipop])
        srcs_lensing=1;
    }
  }
  //Density: the skewers of the sources and the custom maps (the intensity maps read it before the line of sight)
  *send_dens=par->do_cstm || (par->do_srcs && par->do_skewers);
  //Velocities from the LPT: the RSD of the sources, and the lensing shells of the fast lensing, which add
  //the velocity term of the deflection (RETURN_VEL | RETURN_TID in lensing_get_beam_properties)
  *send_vel=par->lpt_vels && (par->do_srcs || par->do_lensing);
  //Potential: tidal tensor (kappa, lensing, lensing of the sources), ISW, and the velocities if not from the LPT
  *send_npot=need_npot && (par->do_kappa || par->do_isw || par->do_lensing || srcs_lensing ||
                           (par->do_srcs && !par->lpt_vels));
}

static void ring_overlapped(ParamCoLoRe *par,int send_dens,int send_vel,int send_npot)
{
  int i,k,nf=0;
  long size_plane=2*((long)par->n_grid)*(par->n_grid/2+1);
  flouble **slot[5],*orig[5],*cur[5],*alt[5];
  long size[5],nbatch[5],nbatch_max=0;
  long off_left=0,off_right=0;
  flouble *hidden_npot=par->grid_npot,*hidden_left=par->slice_left,*hidden_right=par->slice_right;
  flouble *hidden_dens=par->grid_dens,*hidden_velx=par->grid_velx,*hidden_vely=par->grid_vely,*hidden_velz=par->grid_velz;
  if(send_npot) {
    slot[nf]=&(par->grid_npot);
    size[nf]=(par->nz_max+2)*size_plane;
    off_left=par->slice_left-par->grid_npot;
    off_right=par->slice_right-par->grid_npot;
    nf++;
  }
  if(send_dens) {
    slot[nf]=&(par->grid_dens); size[nf]=par->nz_max*size_plane; nf++;
  }
  if(send_vel) {
    slot[nf]=&(par->grid_velx); size[nf]=par->nz_max*size_plane; nf++;
    slot[nf]=&(par->grid_vely); size[nf]=par->nz_max*size_plane; nf++;
    slot[nf]=&(par->grid_velz); size[nf]=par->nz_max*size_plane; nf++;
  }
  for(k=0;k<nf;k++) {
    orig[k]=*(slot[k]);
    cur[k]=orig[k];
    alt[k]=my_malloc(size[k]*sizeof(flouble));
    nbatch[k]=(size[k]+SENDRECV_BATCH-1)/SENDRECV_BATCH;
    if(nbatch[k]>nbatch_max) nbatch_max=nbatch[k];
  }
  MPI_Request *req=my_malloc((2*nf*nbatch_max+1)*sizeof(MPI_Request));
  //The fields that are not sent are not to be read
  if(!send_npot) {
    par->grid_npot=NULL;
    par->slice_left=NULL;
    par->slice_right=NULL;
  }
  if(!send_dens) par->grid_dens=NULL;
  if(!send_vel) {
    par->grid_velx=NULL;
    par->grid_vely=NULL;
    par->grid_velz=NULL;
  }
  for(i=0;i<NNodes;i++) {
    int node_i_am_now=(NodeThis-i+NNodes)%NNodes;
    int nreq=0;
    double t_proc0;
    for(k=0;k<nf;k++)
      *(slot[k])=cur[k];
    if(send_npot) {
      par->slice_left=par->grid_npot+off_left;
      par->slice_right=par->grid_npot+off_right;
    }
    par->nz_here=par->nz_all[node_i_am_now];
    par->iz0_here=par->iz0_all[node_i_am_now];
    if(i<NNodes-1) {
      int b;
      for(k=0;k<nf;k++) {
        for(b=0;b<nbatch[k];b++) {
          long lo=b*(long)SENDRECV_BATCH,n=size[k]-lo;
          if(n>SENDRECV_BATCH) n=SENDRECV_BATCH;
          MPI_Isend(&(cur[k][lo]),(int)n,FLOUBLE_MPI,NodeRight,(int)(k*nbatch_max+b),MPI_COMM_WORLD,&(req[nreq++]));
          MPI_Irecv(&(alt[k][lo]),(int)n,FLOUBLE_MPI,NodeLeft ,(int)(k*nbatch_max+b),MPI_COMM_WORLD,&(req[nreq++]));
        }
      }
    }
    t_proc0=MPI_Wtime();
    beam_operations_one_slab(par);
    if(getenv("COLORE_TIMING")) {
      double t_proc=MPI_Wtime()-t_proc0,t_max,t_sum;
      MPI_Allreduce(&t_proc,&t_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      MPI_Allreduce(&t_proc,&t_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      print_info("  - slab %d processed (processing time over tasks: max %.2f s, mean %.2f s)\n",i,t_max,t_sum/NNodes);
    }
    if(i<NNodes-1) {
      double t_wait0=MPI_Wtime();
      MPI_Waitall(nreq,req,MPI_STATUSES_IGNORE);
      t_ring_mpi+=MPI_Wtime()-t_wait0; //time waiting for the transfer after the processing
      for(k=0;k<nf;k++) {
        flouble *tmp=cur[k];
        cur[k]=alt[k];
        alt[k]=tmp;
      }
    }
  }
  //Back to the original buffers (the content of the fields is not needed after the line of sight)
  for(k=0;k<nf;k++)
    *(slot[k])=orig[k];
  if(send_npot) {
    par->slice_left=par->grid_npot+off_left;
    par->slice_right=par->grid_npot+off_right;
  }
  else {
    par->grid_npot=hidden_npot;
    par->slice_left=hidden_left;
    par->slice_right=hidden_right;
  }
  if(!send_dens) par->grid_dens=hidden_dens;
  if(!send_vel) {
    par->grid_velx=hidden_velx;
    par->grid_vely=hidden_vely;
    par->grid_velz=hidden_velz;
  }
  par->nz_here=par->nz_all[NodeThis];
  par->iz0_here=par->iz0_all[NodeThis];
  for(k=0;k<nf;k++) {
    //whichever of the two buffers is not an original one was allocated here
    if(cur[k]!=orig[k]) free(cur[k]);
    else free(alt[k]);
  }
  free(req);
  if(getenv("COLORE_TIMING"))
    print_info("  - overlapped ring: task 0 waited %.2f s for the transfers after processing\n",t_ring_mpi);
}
#endif //_HAVE_MPI

void get_beam_properties(ParamCoLoRe *par)
{
  print_info("*** Getting LOS information\n");
  if(!par->need_beaming) {
    print_info("  No need!\n\n");
    return;
  }

  if(par->do_kappa)
    kappa_beams_preproc(par);
#ifdef _USE_FAST_LENSING
  if(par->do_lensing)
    lensing_beams_preproc(par);
#endif //_USE_FAST_LENSING
  if(par->do_isw)
    isw_beams_preproc(par);
  if(par->do_srcs)
    srcs_beams_preproc(par);
  if(par->do_imap)
    imap_beams_preproc(par);
  if(par->do_cstm)
    cstm_beams_preproc(par);

  if(NodeThis==0) timer(0);

  int i;
#ifdef _HAVE_MPI
  long size_slice_npot=(par->nz_max+2)*((long)(2*par->n_grid*(par->n_grid/2+1)));
  long size_slice_dens=par->nz_max*((long)(2*par->n_grid*(par->n_grid/2+1)));
  long size_slice_velx=par->nz_max*((long)(2*par->n_grid*(par->n_grid/2+1)));
  long size_slice_vely=par->nz_max*((long)(2*par->n_grid*(par->n_grid/2+1)));
  long size_slice_velz=par->nz_max*((long)(2*par->n_grid*(par->n_grid/2+1)));
  flouble *buffer_sr=my_malloc(size_slice_npot*sizeof(flouble));
  //Without a use for the Newtonian potential do not send it around the ring (one of the five fields
  //exchanged); it is then not even allocated (see allocate_fftw).
  int need_npot=need_newtonian_potential(par);
#endif //_HAVE_MPI

  int use_overlap=0;
#ifdef _HAVE_MPI
  use_overlap=(NNodes>1);
  if(use_overlap) {
    int send_dens,send_vel,send_npot;
    ring_fields_needed(par,need_npot,&send_dens,&send_vel,&send_npot);
    print_info("  Fields sent around the ring: %s%s%s\n",send_dens ? "density " : "",send_vel ? "velocities " : "",
               send_npot ? "potential" : "");
    ring_overlapped(par,send_dens,send_vel,send_npot);
  }
#endif //_HAVE_MPI

  for(i=0;i<NNodes && !use_overlap;i++) {
    int node_i_am_now=(NodeThis-i-1+NNodes)%NNodes;
#ifdef _HAVE_MPI
#ifdef _DEBUG
    print_info("Communication %d, MPI task %d is now MPI task %d\n",i,NodeThis,node_i_am_now);
#endif //_DEBUG
    if(need_npot)
      mpi_sendrecv_wrap(par->grid_npot,buffer_sr,size_slice_npot,i);
    mpi_sendrecv_wrap(par->grid_dens,buffer_sr,size_slice_dens,i);
    if (par->lpt_vels) {
        mpi_sendrecv_wrap(par->grid_velx,buffer_sr,size_slice_velx,i);
        mpi_sendrecv_wrap(par->grid_vely,buffer_sr,size_slice_vely,i);
        mpi_sendrecv_wrap(par->grid_velz,buffer_sr,size_slice_velz,i);
        }
#endif //_HAVE_MPI
    par->nz_here=par->nz_all[node_i_am_now];
    par->iz0_here=par->iz0_all[node_i_am_now];
    if(getenv("COLORE_TIMING")) print_info("  - slab %d exchanged (task 0, total so far: MPI %.2f s, copy %.2f s)\n",i,t_ring_mpi,t_ring_copy);
#ifdef _HAVE_MPI
    double t_proc0=MPI_Wtime();
#endif //_HAVE_MPI

    beam_operations_one_slab(par);
    if(getenv("COLORE_TIMING")) {
#ifdef _HAVE_MPI
      double t_proc=MPI_Wtime()-t_proc0,t_max,t_sum;
      MPI_Allreduce(&t_proc,&t_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      MPI_Allreduce(&t_proc,&t_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      print_info("  - slab %d processed (processing time over tasks: max %.2f s, mean %.2f s)\n",i,t_max,t_sum/NNodes);
#else //_HAVE_MPI
      print_info("  - slab %d processed\n",i);
#endif //_HAVE_MPI
    }
  }
#ifdef _HAVE_MPI
  free(buffer_sr);
#endif //_HAVE_MPI

  if(par->do_kappa)
    kappa_beams_postproc(par);
#ifdef _USE_FAST_LENSING
  if(par->do_lensing)
    lensing_beams_postproc(par); 
#endif //_USE_FAST_LENSING
  if(par->do_isw)
    isw_beams_postproc(par);
  if(par->do_srcs)
    srcs_beams_postproc(par);
  if(par->do_imap)
    imap_beams_postproc(par);
  if(par->do_cstm)
    cstm_beams_postproc(par);

  if(NodeThis==0) timer(2);
  print_info("\n");
}
