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
			flouble *d,flouble v[3],flouble t[6],flouble *pdot,
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

  //Get velocity
  if(flag_return & RETURN_VEL) {
    v[0]=par->grid_npot[ix_hi+iy_0+iz_0]-par->grid_npot[ix_lo+iy_0+iz_0];
    v[1]=par->grid_npot[ix_0+iy_hi+iz_0]-par->grid_npot[ix_0+iy_lo+iz_0];
    if(iz==0)
      v[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->slice_left[ix_0+iy_0];
    else if(iz==par->nz_here-1)
      v[2]=par->slice_right[ix_0+iy_0]-par->grid_npot[ix_0+iy_0+iz_lo];
    else
      v[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->grid_npot[ix_0+iy_0+iz_lo];
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
			  flouble *d,flouble v[3],flouble t[6],flouble *pd,
			  int flag_return,int interp_type)
{
  long ix0[3];
  double h0x[3];
  int ax,added_anything=0;

  //Initialize output
  if(flag_return & RETURN_DENS)
    *d=0;
  
  if(flag_return & RETURN_VEL) {
    for(ax=0;ax<3;ax++)
      v[ax]=0;
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
      flouble d_000,v_000[3],t_000[6],pd_000;
      flouble w_000=h0x[2]*h0x[1]*h0x[0];
      
      added_anything=1;
      get_element(par,ix0[0],ix0[1],ix0[2],&d_000,v_000,t_000,&pd_000,flag_return);
      if(flag_return & RETURN_DENS)
	*d+=d_000*w_000;
      if(flag_return & RETURN_VEL) {
	for(ax=0;ax<3;ax++)
	  v[ax]+=v_000[ax]*w_000;
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
      flouble d_000,v_000[3],t_000[6],pd_000;
      flouble d_001,v_001[3],t_001[6],pd_001;
      flouble d_010,v_010[3],t_010[6],pd_010;
      flouble d_011,v_011[3],t_011[6],pd_011;
      flouble w_000=h1x[2]*h1x[1]*h1x[0];
      flouble w_001=h1x[2]*h1x[1]*h0x[0];
      flouble w_010=h1x[2]*h0x[1]*h1x[0];
      flouble w_011=h1x[2]*h0x[1]*h0x[0];
      
      added_anything=1;
      get_element(par,ix0[0],ix0[1],ix0[2],&d_000,v_000,t_000,&pd_000,flag_return);
      get_element(par,ix1[0],ix0[1],ix0[2],&d_001,v_001,t_001,&pd_001,flag_return);
      get_element(par,ix0[0],ix1[1],ix0[2],&d_010,v_010,t_010,&pd_010,flag_return);
      get_element(par,ix1[0],ix1[1],ix0[2],&d_011,v_011,t_011,&pd_011,flag_return);
      if(flag_return & RETURN_DENS)
	*d+=(d_000*w_000+d_001*w_001+d_010*w_010+d_011*w_011);
      if(flag_return & RETURN_VEL) {
	for(ax=0;ax<3;ax++)
	  v[ax]+=(v_000[ax]*w_000+v_001[ax]*w_001+v_010[ax]*w_010+v_011[ax]*w_011);
      }
      if(flag_return & RETURN_TID) {
	for(ax=0;ax<6;ax++)
	  t[ax]+=(t_000[ax]*w_000+t_001[ax]*w_001+t_010[ax]*w_010+t_011[ax]*w_011);
      }
      if(flag_return & RETURN_PDOT)
	*pd+=(pd_000*w_000+pd_001*w_001+pd_010*w_010+pd_011*w_011);
    }
    if((ix1[2]>=0) && (ix1[2]<par->nz_here)) {
      flouble d_100,v_100[3],t_100[6],pd_100;
      flouble d_101,v_101[3],t_101[6],pd_101;
      flouble d_110,v_110[3],t_110[6],pd_110;
      flouble d_111,v_111[3],t_111[6],pd_111;
      flouble w_100=h0x[2]*h1x[1]*h1x[0];
      flouble w_101=h0x[2]*h1x[1]*h0x[0];
      flouble w_110=h0x[2]*h0x[1]*h1x[0];
      flouble w_111=h0x[2]*h0x[1]*h0x[0];
      
      added_anything=1;
      get_element(par,ix0[0],ix0[1],ix1[2],&d_100,v_100,t_100,&pd_100,flag_return);
      get_element(par,ix1[0],ix0[1],ix1[2],&d_101,v_101,t_101,&pd_101,flag_return);
      get_element(par,ix0[0],ix1[1],ix1[2],&d_110,v_110,t_110,&pd_110,flag_return);
      get_element(par,ix1[0],ix1[1],ix1[2],&d_111,v_111,t_111,&pd_111,flag_return);
      if(flag_return & RETURN_DENS)
	*d+=(d_100*w_100+d_101*w_101+d_110*w_110+d_111*w_111);
      if(flag_return & RETURN_VEL) {
	for(ax=0;ax<3;ax++)
	  v[ax]+=(v_100[ax]*w_100+v_101[ax]*w_101+v_110[ax]*w_110+v_111[ax]*w_111);
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

#ifdef _HAVE_MPI 
static void mpi_sendrecv_wrap(flouble *data,flouble *buff,long count,int tag)
{
  // still need to compile even if never called
#define SENDRECV_BATCH 1073741824
  int remainder;
  long i_sofar=0;
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
  memcpy(data,buff,count*sizeof(flouble));
}
#endif //_HAVE_MPI

void get_beam_properties(ParamCoLoRe *par)
{
  print_info("*** Getting LOS information\n");
  if(!par->need_beaming) {
    print_info("  No need!\n\n");
    return;
  }

  if(par->do_srcs)
    srcs_beams_preproc(par);
  if(par->do_imap)
    imap_beams_preproc(par);
  if(par->do_kappa)
    kappa_beams_preproc(par);
  if(par->do_isw)
    isw_beams_preproc(par);
  
  if(NodeThis==0) timer(0);

  int i;
#ifdef _HAVE_MPI
  long size_slice_npot=(par->nz_max+2)*((long)(2*par->n_grid*(par->n_grid/2+1)));
  long size_slice_dens=par->nz_max*((long)(2*par->n_grid*(par->n_grid/2+1)));
  flouble *buffer_sr=my_malloc(size_slice_npot*sizeof(flouble));
#endif //_HAVE_MPI

  for(i=0;i<NNodes;i++) {
    int node_i_am_now=(NodeThis-i-1+NNodes)%NNodes;
#ifdef _HAVE_MPI
#ifdef _DEBUG
    print_info("Communication %d, Node %d is now Node %d\n",i,NodeThis,node_i_am_now);
#endif //_DEBUG
    mpi_sendrecv_wrap(par->grid_npot,buffer_sr,size_slice_npot,i);
    mpi_sendrecv_wrap(par->grid_dens,buffer_sr,size_slice_dens,i);

#endif //_HAVE_MPI
    par->nz_here=par->nz_all[node_i_am_now];
    par->iz0_here=par->iz0_all[node_i_am_now];

    if(par->do_srcs)
      srcs_get_beam_properties(par);
    if(par->do_imap)
      imap_get_beam_properties(par);
    if(par->do_kappa)
      kappa_get_beam_properties(par);
    if(par->do_isw)
      isw_get_beam_properties(par);
  }
#ifdef _HAVE_MPI
  free(buffer_sr);
#endif //_HAVE_MPI

  if(par->do_srcs)
    srcs_beams_postproc(par);
  if(par->do_imap)
    imap_beams_postproc(par);
  if(par->do_kappa)
    kappa_beams_postproc(par);
  if(par->do_isw)
    isw_beams_postproc(par);
  
  if(NodeThis==0) timer(2);
  print_info("\n");
}
