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

#define INTERP_NGP 0
#define INTERP_CIC 1
#ifndef INTERP_TYPE
#define INTERP_TYPE INTERP_NGP
#endif //INTERP_TYPE

#define IND_XX 0
#define IND_XY 1
#define IND_XZ 2
#define IND_YY 3
#define IND_YZ 4
#define IND_ZZ 5

static inline double get_cosine(double index,double dx)
{
#if PIXTYPE==PT_CEA
  return index*dx-1;
#elif PIXTYPE==PT_CAR
  return cos(M_PI-index*dx);
#endif //PIXTYPE
}

static void get_element(ParamCoLoRe *par,int ix,int iy,int iz,
			flouble *d,flouble v[3],flouble t[6])
{
  int ngx=2*(par->n_grid/2+1);
  lint iz_hi=iz+1,iz_lo=iz-1,iz_0=iz;
  lint iy_hi=iy+1,iy_lo=iy-1,iy_0=iy;
  lint ix_hi=ix+1,ix_lo=ix-1,ix_0=ix;
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
  *d=par->grid_dens[ix_0+iy_0+iz_0];

  //Get velocity
  v[0]=par->grid_npot[ix_hi+iy_0+iz_0]-par->grid_npot[ix_lo+iy_0+iz_0];
  v[1]=par->grid_npot[ix_0+iy_hi+iz_0]-par->grid_npot[ix_0+iy_lo+iz_0];
  if(iz==0)
    v[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->slice_left[ix_0+iy_0];
  else if(iz==par->nz_here-1)
    v[2]=par->slice_right[ix_0+iy_0]-par->grid_npot[ix_0+iy_0+iz_lo];
  else
    v[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->grid_npot[ix_0+iy_0+iz_lo];

  //Get tidal tensor
  if(par->do_lensing) {
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

void pixelize(ParamCoLoRe *par)
{
  print_info("*** Pixelizing cartesian grids\n");
  if(NodeThis==0) timer(0);

  int i;
  lint size_slice_npot=(par->nz_max+2)*((lint)(2*par->n_grid*(par->n_grid/2+1)));
  lint size_slice_dens=par->nz_max*((lint)(2*par->n_grid*(par->n_grid/2+1)));

  for(i=0;i<NNodes;i++) {
    int ib;
    int node_i_am_now=(NodeThis-i-1+NNodes)%NNodes;
#ifdef _HAVE_MPI
#ifdef _DEBUG
    print_info("Communication %d, Node %d is now Node %d\n",i,NodeThis,node_i_am_now);
#endif //_DEBUG
    MPI_Sendrecv_replace(par->grid_npot,size_slice_npot,FLOUBLE_MPI,
			 NodeRight,i,NodeLeft,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Sendrecv_replace(par->grid_dens,size_slice_dens,FLOUBLE_MPI,
			 NodeRight,i,NodeLeft,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
#endif //_HAVE_MPI
    par->nz_here=par->nz_all[node_i_am_now];
    par->iz0_here=par->iz0_all[node_i_am_now];

    for(ib=0;ib<par->n_beams_here;ib++) {
      OnionInfo *oi=par->oi_beams[ib];
#ifdef _HAVE_OMP
#pragma omp parallel default(none) shared(par,oi,ib)
#endif //_HAVE_OMP
      {
	int ir;
	flouble idx=par->n_grid/par->l_box;
	double factor_vel=-par->fgrowth_0/(1.5*par->hubble_0*par->OmegaM);

#ifdef _HAVE_OMP
#pragma omp for nowait schedule(dynamic)
#endif //_HAVE_OMP
	for(ir=0;ir<oi->nr;ir++) {
	  int icth;
	  double dr=oi->rf_arr[ir]-oi->r0_arr[ir];
#if PIXTYPE==PT_CEA
	  flouble dcth=2.0/oi->nside_arr[ir];
	  flouble dcth_sub=dcth/FAC_CART2SPH_NSUB;
#elif PIXTYPE==PT_CAR
	  flouble dth=M_PI/oi->nside_arr[ir];
#endif //PIXTYPE
	  flouble dphi=M_PI/oi->nside_arr[ir];
	  flouble dr_sub=dr/FAC_CART2SPH_NSUB;
	  flouble dphi_sub=dphi/FAC_CART2SPH_NSUB;
	  int ncth=oi->icthf_arr[ir]-oi->icth0_arr[ir]+1;
	  int nphi=oi->iphif_arr[ir]-oi->iphi0_arr[ir]+1;
	  flouble r0=oi->r0_arr[ir];
	  //	  printf("%d %d %d %d\n",ib,ir,ncth,nphi);

	  for(icth=0;icth<ncth;icth++) {
	    int iphi;
	    int index_cth=nphi*icth;
#if PIXTYPE==PT_CEA
	    flouble cth0=get_cosine(oi->icth0_arr[ir]+icth+0.0,dcth);
#elif PIXTYPE==PT_CAR
	    flouble cth0=get_cosine(oi->icth0_arr[ir]+icth+0.0,dth);
	    flouble dcth=get_cosine(oi->icth0_arr[ir]+icth+1.0,dth)-cth0;
	    flouble dcth_sub=dcth/FAC_CART2SPH_NSUB;
#endif //PIXTYPE
	    flouble cth_h=cth0+dcth*0.5;
	    flouble sth_h=sqrt(1-cth_h*cth_h);
	    for(iphi=0;iphi<nphi;iphi++) {
	      int ir2;
	      flouble u[3];
	      int added_anything=0;
	      int index=index_cth+iphi;
	      flouble phi0=(oi->iphi0_arr[ir]+iphi)*dphi;
	      flouble phi_h=phi0+dphi*0.5;
	      flouble cph_h=cos(phi_h),sph_h=sin(phi_h);
	      flouble d=0,v[3]={0,0,0},t[6]={0,0,0,0,0,0};
	      u[0]=sth_h*cph_h; u[1]=sth_h*sph_h; u[2]=cth_h;

	      //Make sub-voxels
	      for(ir2=0;ir2<FAC_CART2SPH_NSUB;ir2++) {
		int icth2;
		flouble r=r0+(ir2+0.5)*dr_sub;
		for(icth2=0;icth2<FAC_CART2SPH_NSUB;icth2++) {
		  int iphi2;
		  flouble cth=cth0+(icth2+0.5)*dcth_sub;
		  flouble sth=sqrt(1-cth*cth);
		  for(iphi2=0;iphi2<FAC_CART2SPH_NSUB;iphi2++) {
		    flouble phi=phi0+(iphi2+0.5)*dphi_sub;
		    int ax;
		    lint ix0[3];
		    flouble x[3],h0x[3];
		    x[0]=r*sth*cos(phi); x[1]=r*sth*sin(phi); x[2]=r*cth;
#if INTERP_TYPE==INTERP_NGP
		    //Trilinear interpolation
		    for(ax=0;ax<3;ax++) {
		      x[ax]+=par->pos_obs[ax];
		      ix0[ax]=(int)(x[ax]*idx+0.5);
		      if(ix0[ax]>=par->n_grid)
			ix0[ax]-=par->n_grid;
		      else if(ix0[ax]<0)
			ix0[ax]+=par->n_grid;
		      h0x[ax]=1./FAC_CART2SPH_NSUB;
		    }
		    ix0[2]-=par->iz0_here;
		    
		    if((ix0[2]>=0) && (ix0[2]<par->nz_here)) {
		      flouble d_000,v_000[3],t_000[6];
		      flouble w_000=h0x[2]*h0x[1]*h0x[0];
		      
		      added_anything=1;
		      get_element(par,ix0[0],ix0[1],ix0[2],&d_000,v_000,t_000);
		      d+=d_000*w_000;
		      for(ax=0;ax<3;ax++)
			v[ax]+=v_000[ax]*w_000;
		      if(par->do_lensing) {
			for(ax=0;ax<6;ax++)
			  t[ax]+=t_000[ax]*w_000;
		      }
		    }
#elif INTERP_TYPE==INTERP_CIC
		    lint ix1[3];
		    flouble h1x[3];

		    //Trilinear interpolation
		    for(ax=0;ax<3;ax++) {
		      x[ax]+=par->pos_obs[ax];
		      ix0[ax]=(int)(x[ax]*idx);
		      h0x[ax]=x[ax]*idx-ix0[ax];
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
		      h0x[ax]/=FAC_CART2SPH_NSUB;
		      h1x[ax]/=FAC_CART2SPH_NSUB;
		    }
		    ix0[2]-=par->iz0_here;
		    ix1[2]-=par->iz0_here;
		    
		    if((ix0[2]>=0) && (ix0[2]<par->nz_here)) {
		      flouble d_000,v_000[3],t_000[6];
		      flouble d_001,v_001[3],t_001[6];
		      flouble d_010,v_010[3],t_010[6];
		      flouble d_011,v_011[3],t_011[6];
		      flouble w_000=h1x[2]*h1x[1]*h1x[0];
		      flouble w_001=h1x[2]*h1x[1]*h0x[0];
		      flouble w_010=h1x[2]*h0x[1]*h1x[0];
		      flouble w_011=h1x[2]*h0x[1]*h0x[0];
		      
		      added_anything=1;
		      get_element(par,ix0[0],ix0[1],ix0[2],&d_000,v_000,t_000);
		      get_element(par,ix1[0],ix0[1],ix0[2],&d_001,v_001,t_001);
		      get_element(par,ix0[0],ix1[1],ix0[2],&d_010,v_010,t_010);
		      get_element(par,ix1[0],ix1[1],ix0[2],&d_011,v_011,t_011);
		      d+=(d_000*w_000+d_001*w_001+d_010*w_010+d_011*w_011);
		      for(ax=0;ax<3;ax++)
			v[ax]+=(v_000[ax]*w_000+v_001[ax]*w_001+v_010[ax]*w_010+v_011[ax]*w_011);
		      if(par->do_lensing) {
			for(ax=0;ax<6;ax++)
			  t[ax]+=(t_000[ax]*w_000+t_001[ax]*w_001+t_010[ax]*w_010+t_011[ax]*w_011);
		      }
		    }
		    if((ix1[2]>=0) && (ix1[2]<par->nz_here)) {
		      flouble d_100,v_100[3],t_100[6];
		      flouble d_101,v_101[3],t_101[6];
		      flouble d_110,v_110[3],t_110[6];
		      flouble d_111,v_111[3],t_111[6];
		      flouble w_100=h0x[2]*h1x[1]*h1x[0];
		      flouble w_101=h0x[2]*h1x[1]*h0x[0];
		      flouble w_110=h0x[2]*h0x[1]*h1x[0];
		      flouble w_111=h0x[2]*h0x[1]*h0x[0];
		      
		      added_anything=1;
		      get_element(par,ix0[0],ix0[1],ix1[2],&d_100,v_100,t_100);
		      get_element(par,ix1[0],ix0[1],ix1[2],&d_101,v_101,t_101);
		      get_element(par,ix0[0],ix1[1],ix1[2],&d_110,v_110,t_110);
		      get_element(par,ix1[0],ix1[1],ix1[2],&d_111,v_111,t_111);
		      d+=(d_100*w_100+d_101*w_101+d_110*w_110+d_111*w_111);
		      for(ax=0;ax<3;ax++)
			v[ax]+=(v_100[ax]*w_100+v_101[ax]*w_101+v_110[ax]*w_110+v_111[ax]*w_111);
		      if(par->do_lensing) {
			for(ax=0;ax<6;ax++)
			  t[ax]+=(t_100[ax]*w_100+t_101[ax]*w_101+t_110[ax]*w_110+t_111[ax]*w_111);
		      }
		    }
#endif //INTERP_TYPE
		  }
		}
	      }

	      if(added_anything) {
		par->dens_beams[ib][ir][index]+=d;
		par->vrad_beams[ib][ir][index]+=factor_vel*0.5*idx*(v[0]*u[0]+v[1]*u[1]+v[2]*u[2]);
		if(par->do_lensing) {
		  par->p_xx_beams[ib][ir][index]+=idx*idx*
		    (cth_h*cth_h*(t[IND_XX]*cph_h*cph_h+2*t[IND_XY]*cph_h*sph_h+t[IND_YY]*sph_h*sph_h)+
		     t[IND_ZZ]*sth_h*sth_h-2*cth_h*sth_h*(t[IND_XZ]*cph_h+t[IND_YZ]*sph_h));
		  par->p_xy_beams[ib][ir][index]+=idx*idx*
		    (t[IND_XY]*(cph_h*cph_h-sph_h*sph_h)*cth_h+t[IND_XZ]*sph_h*sth_h-
		     cph_h*((t[IND_XX]-t[IND_YY])*cth_h*sph_h+t[IND_YZ]*sth_h));
		  par->p_yy_beams[ib][ir][index]+=idx*idx*
		    (t[IND_XX]*sph_h*sph_h+t[IND_YY]*cph_h*cph_h-2*t[IND_XY]*cph_h*sph_h);
		}
	      }
	    }
	  }
	} //end omp for
      } //end omp parallel
    }
  }

  if(NodeThis==0) timer(2);
  end_fftw(par);
  print_info("\n");
}
