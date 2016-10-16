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

#define CPY_TASK_DENS 0
#define CPY_TASK_VRAD 1
#define CPY_TASK_P_XX 2
#define CPY_TASK_P_XY 3
#define CPY_TASK_P_YY 4

static inline double get_cosine(double index,double dx)
{
#if PIXTYPE==PT_CEA
  return index*dx-1;
#elif PIXTYPE==PT_CAR
  return cos(M_PI-index*dx);
#endif //PIXTYPE
}

//////
// Copies new observable quantity into dummy array
static void copy_to_dum(ParamCoLoRe *par,int cpy_task)
{
#ifdef _DEBUG
  print_info("Copying to dum\n");
  if(NodeThis==0) timer(0);
#endif //_DEBUG
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,cpy_task)
#endif //_HAVE_OMP
  {
    double dx=par->l_box/par->n_grid;
    double idx=1./dx;
    int iz;
    int ngx=2*(par->n_grid/2+1);
    double factor_vel=-par->fgrowth_0/(1.5*par->hubble_0*par->OmegaM);
    gsl_matrix *mat_tij=NULL,*mat_basis=NULL,*mat_dum=NULL;
    if((cpy_task==CPY_TASK_P_XX) || (cpy_task==CPY_TASK_P_XY) || (cpy_task==CPY_TASK_P_YY)) {
      mat_tij=gsl_matrix_alloc(3,3);
      mat_basis=gsl_matrix_alloc(3,3);
      mat_dum=gsl_matrix_alloc(3,3);
    }

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint iz_hi=iz+1;
      lint iz_lo=iz-1;
      lint iz_0=iz;
      double z=dx*(iz+par->iz0_here+0.5)-par->pos_obs[2];
      if(iz==0) iz_lo=par->n_grid-1;
      if(iz==par->n_grid-1) iz_hi=0;
      iz_hi*=ngx*par->n_grid;
      iz_lo*=ngx*par->n_grid;
      iz_0*=ngx*par->n_grid;
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint iy_hi=iy+1;
	lint iy_lo=iy-1;
	lint iy_0=iy;
	double y=dx*(iy+0.5)-par->pos_obs[1];
	if(iy==0) iy_lo=par->n_grid-1;
	if(iy==par->n_grid-1) iy_hi=0;
	iy_hi*=ngx;
	iy_lo*=ngx;
	iy_0*=ngx;
	for(ix=0;ix<par->n_grid;ix++) {
	  lint ix_hi=ix+1;
	  lint ix_lo=ix-1;
	  lint ix_0=ix;
	  double x=dx*(ix+0.5)-par->pos_obs[0];
	  if(ix==0) ix_lo=par->n_grid-1;
	  if(ix==par->n_grid-1) ix_hi=0;
	  if(cpy_task==CPY_TASK_DENS) {
	    par->grid_dumm[ix_0+iy_0+iz_0]=par->grid_dens[ix_0+iy_0+iz_0];
	  }
	  else if(cpy_task==CPY_TASK_VRAD) {
	    double vel[3];
	    double irr=1./sqrt(x*x+y*y+z*z);

	    vel[0]=par->grid_npot[ix_hi+iy_0+iz_0]-par->grid_npot[ix_lo+iy_0+iz_0];
	    vel[1]=par->grid_npot[ix_0+iy_hi+iz_0]-par->grid_npot[ix_0+iy_lo+iz_0];
	    if(iz==0)
	      vel[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->slice_left[ix_0+iy_0];
	    else if(iz==par->nz_here-1)
	      vel[2]=par->slice_right[ix_0+iy_0]-par->grid_npot[ix_0+iy_0+iz_lo];
	    else
	      vel[2]=par->grid_npot[ix_0+iy_0+iz_hi]-par->grid_npot[ix_0+iy_0+iz_lo];

	    par->grid_dumm[ix_0+iy_0+iz_0]=factor_vel*0.5*idx*irr*(vel[0]*x+vel[1]*y+vel[2]*z);
	  }
	  else if((cpy_task==CPY_TASK_P_XX) || (cpy_task==CPY_TASK_P_XY) || (cpy_task==CPY_TASK_P_YY)) {
	    double cth,sth,cph,sph;
	    double r=sqrt(x*x+y*y+z*z);
	    double txx,txy,txz,tyy,tyz,tzz;

	    //Set rotation matrix
	    if(r==0) {
	      cth=1; cph=1; sph=0; sth=0;
	    }
	    else {
	      cth=z/r;
	      sth=sqrt(1-cth*cth);
	      if(sth==0) {
		cph=1;
		sph=0;
	      }
	      else {
		cph=x/(r*sth);
		sph=y/(r*sth);
	      }
	    }
	    gsl_matrix_set(mat_basis,0,0,sth*cph);
	    gsl_matrix_set(mat_basis,0,1,sth*sph);
	    gsl_matrix_set(mat_basis,0,2,cth);
	    gsl_matrix_set(mat_basis,1,0,cth*cph);
	    gsl_matrix_set(mat_basis,1,1,cth*sph);
	    gsl_matrix_set(mat_basis,1,2,-sth);
	    gsl_matrix_set(mat_basis,2,0,-sph);
	    gsl_matrix_set(mat_basis,2,1,cph);
	    gsl_matrix_set(mat_basis,2,2,0);

	    //Set tidal tensor
	    txx=(par->grid_npot[ix_hi+iy_0+iz_0]+par->grid_npot[ix_lo+iy_0+iz_0]-2*par->grid_npot[ix_0+iy_0+iz_0]);
	    tyy=(par->grid_npot[ix_0+iy_hi+iz_0]+par->grid_npot[ix_0+iy_lo+iz_0]-2*par->grid_npot[ix_0+iy_0+iz_0]);
	    txy=0.25*(par->grid_npot[ix_hi+iy_hi+iz_0]+par->grid_npot[ix_lo+iy_lo+iz_0]-
		      par->grid_npot[ix_hi+iy_lo+iz_0]-par->grid_npot[ix_lo+iy_hi+iz_0]);
	    if(iz==0) {
	      tzz=(par->grid_npot[ix_0+iy_0+iz_hi]+par->slice_left[ix_0+iy_0]-2*par->grid_npot[ix_0+iy_0+iz_0]);
	      txz=0.25*(par->grid_npot[ix_hi+iy_0+iz_hi]+par->slice_left[ix_lo+iy_0]-
			par->slice_left[ix_hi+iy_0]-par->grid_npot[ix_lo+iy_0+iz_hi]);
	      tyz=0.25*(par->grid_npot[ix_0+iy_hi+iz_hi]+par->slice_left[ix_0+iy_lo]-
			par->slice_left[ix_0+iy_hi]-par->grid_npot[ix_0+iy_lo+iz_hi]);
	    }
	    else if(iz==par->nz_here-1) {
	      tzz=(par->slice_right[ix_0+iy_0]+par->grid_npot[ix_0+iy_0+iz_lo]-2*par->grid_npot[ix_0+iy_0+iz_0]);
	      txz=0.25*(par->slice_right[ix_hi+iy_0]+par->grid_npot[ix_lo+iy_0+iz_lo]-
			par->grid_npot[ix_hi+iy_0+iz_lo]-par->slice_right[ix_lo+iy_0]);
	      tyz=0.25*(par->slice_right[ix_0+iy_hi]+par->grid_npot[ix_0+iy_lo+iz_lo]-
			par->grid_npot[ix_0+iy_hi+iz_lo]-par->slice_right[ix_0+iy_lo]);
	    }
	    else {
	      tzz=(par->grid_npot[ix_0+iy_0+iz_hi]+par->grid_npot[ix_0+iy_0+iz_lo]-2*par->grid_npot[ix_0+iy_0+iz_0]);
	      txz=0.25*(par->grid_npot[ix_hi+iy_0+iz_hi]+par->grid_npot[ix_lo+iy_0+iz_lo]-
			par->grid_npot[ix_hi+iy_0+iz_lo]-par->grid_npot[ix_lo+iy_0+iz_hi]);
	      tyz=0.25*(par->grid_npot[ix_0+iy_hi+iz_hi]+par->grid_npot[ix_0+iy_lo+iz_lo]-
			par->grid_npot[ix_0+iy_hi+iz_lo]-par->grid_npot[ix_0+iy_lo+iz_hi]);
	    }
	    gsl_matrix_set(mat_tij,0,0,txx);
	    gsl_matrix_set(mat_tij,0,1,txy);
	    gsl_matrix_set(mat_tij,0,2,txz);
	    gsl_matrix_set(mat_tij,1,0,txy);
	    gsl_matrix_set(mat_tij,1,1,tyy);
	    gsl_matrix_set(mat_tij,1,2,tyz);
	    gsl_matrix_set(mat_tij,2,0,txz);
	    gsl_matrix_set(mat_tij,2,1,tyz);
	    gsl_matrix_set(mat_tij,2,2,tzz);

	    //Rotate to spherical coordinates
	    gsl_blas_dgemm(CblasNoTrans,CblasTrans  ,1.,mat_tij  ,mat_basis,0,mat_dum);
	    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,mat_basis,mat_dum  ,0,mat_tij);

	    if(cpy_task==CPY_TASK_P_XX)
	      par->grid_dumm[ix_0+iy_0+iz_0]=idx*idx*gsl_matrix_get(mat_tij,1,1);
	    if(cpy_task==CPY_TASK_P_XY)
	      par->grid_dumm[ix_0+iy_0+iz_0]=idx*idx*gsl_matrix_get(mat_tij,1,2);
	    if(cpy_task==CPY_TASK_P_YY)
	      par->grid_dumm[ix_0+iy_0+iz_0]=idx*idx*gsl_matrix_get(mat_tij,2,2);
	  }
	  else
	    report_error(1,"Wrong task %d\n",cpy_task);
	}
      }
    } // end omp for
    if((cpy_task==CPY_TASK_P_XX) || (cpy_task==CPY_TASK_P_XY) || (cpy_task==CPY_TASK_P_YY)) {
      gsl_matrix_free(mat_tij);
      gsl_matrix_free(mat_basis);
      gsl_matrix_free(mat_dum);
    }    
  } // end omp parallel
#ifdef _DEBUG
  if(NodeThis==0) timer(2);
#endif //_DEBUG
}

static void interpolate_to_slice(ParamCoLoRe *par,OnionInfo *oi,flouble *grid,flouble **slices)
{
#ifdef _DEBUG
  print_info("Interpolating to slices\n");
  if(NodeThis==0) timer(0);
#endif //_DEBUG

#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par,oi,grid,slices)
#endif //_HAVE_OMP
  {
    int ir;
    flouble idx=par->n_grid/par->l_box;
    lint ngx=2*(par->n_grid/2+1);

#ifdef _HAVE_OMP
#pragma omp for schedule(dynamic)
#endif //_HAVE_OMP
    for(ir=0;ir<oi->nr;ir++) {
      if(oi->num_pix[ir]>0) {
	int icth;
	flouble dr=oi->rf_arr[ir]-oi->r0_arr[ir];
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
	  for(iphi=0;iphi<nphi;iphi++) {
	    int ir2;
	    flouble phi0=(oi->iphi0_arr[ir]+iphi)*dphi;
	    int index=index_cth+iphi;
	    flouble arr_add=0;

	    //Make sub-voxels
	    for(ir2=0;ir2<FAC_CART2SPH_NSUB;ir2++) {
	      int icth2;
	      flouble r=r0+(ir2+0.5)*dr_sub;
	      for(icth2=0;icth2<FAC_CART2SPH_NSUB;icth2++) {
		int iphi2;
		flouble cth=cth0+(icth2+0.5)*dcth_sub;
		flouble sth=sqrt(1-cth*cth);
		for(iphi2=0;iphi2<FAC_CART2SPH_NSUB;iphi2++) {
		  int ax;
		  lint ix1[3],ix0[3];
		  flouble x[3],h0x[3],h1x[3];
		  flouble phi=phi0+(iphi2+0.5)*dphi_sub;
		  x[0]=r*sth*cos(phi);
		  x[1]=r*sth*sin(phi);
		  x[2]=r*cth;

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
		  }
		  ix0[2]-=par->iz0_here;
		  ix1[2]-=par->iz0_here;
		  if((ix0[2]>=0) && (ix0[2]<par->nz_here)) {
		    arr_add+=
		      grid[ix0[2]*ngx*par->n_grid+ix0[1]*ngx+ix0[0]]*h1x[2]*h1x[1]*h1x[0]+
		      grid[ix0[2]*ngx*par->n_grid+ix0[1]*ngx+ix1[0]]*h1x[2]*h1x[1]*h0x[0]+
		      grid[ix0[2]*ngx*par->n_grid+ix1[1]*ngx+ix0[0]]*h1x[2]*h0x[1]*h1x[0]+
		      grid[ix0[2]*ngx*par->n_grid+ix1[1]*ngx+ix1[0]]*h1x[2]*h0x[1]*h0x[0];
		  }
		  if((ix1[2]>=0) && (ix1[2]<par->nz_here)) {
		    arr_add+=
		      grid[ix1[2]*ngx*par->n_grid+ix0[1]*ngx+ix0[0]]*h0x[2]*h1x[1]*h1x[0]+
		      grid[ix1[2]*ngx*par->n_grid+ix0[1]*ngx+ix1[0]]*h0x[2]*h1x[1]*h0x[0]+
		      grid[ix1[2]*ngx*par->n_grid+ix1[1]*ngx+ix0[0]]*h0x[2]*h0x[1]*h1x[0]+
		      grid[ix1[2]*ngx*par->n_grid+ix1[1]*ngx+ix1[0]]*h0x[2]*h0x[1]*h0x[0];
		  }
		}
	      }
	    }
	    slices[ir][index]+=arr_add/(FAC_CART2SPH_NSUB*FAC_CART2SPH_NSUB*FAC_CART2SPH_NSUB);
	  }
	}
      }
    } //end omp for
  } //end omp parallel
#ifdef _DEBUG
  if(NodeThis==0) timer(2);
#endif //_DEBUG
}

static void compute_slices(ParamCoLoRe *par,int task)
{
  copy_to_dum(par,task);
  if(task==CPY_TASK_DENS)
    interpolate_to_slice(par,par->oi_slices,par->grid_dumm,par->dens_slices);
  else if(task==CPY_TASK_VRAD)
    interpolate_to_slice(par,par->oi_slices,par->grid_dumm,par->vrad_slices);
  else if(task==CPY_TASK_P_XX)
    interpolate_to_slice(par,par->oi_slices,par->grid_dumm,par->p_xx_slices);
  else if(task==CPY_TASK_P_XY)
    interpolate_to_slice(par,par->oi_slices,par->grid_dumm,par->p_xy_slices);
  else if(task==CPY_TASK_P_YY)
    interpolate_to_slice(par,par->oi_slices,par->grid_dumm,par->p_yy_slices);
  else
    report_error(1,"Wrong task %d\n",task);
}

static void communicate_onion_slices(OnionInfo *oi_send,OnionInfo *oi_recv,MPI_Status *stat)
{
#ifdef _HAVE_MPI
  MPI_Sendrecv(oi_send->iphi0_arr,oi_send->nr,MPI_INT,NodeRight,1,
	       oi_recv->iphi0_arr,oi_send->nr,MPI_INT,NodeLeft ,1,
	       MPI_COMM_WORLD,stat);
  MPI_Sendrecv(oi_send->iphif_arr,oi_send->nr,MPI_INT,NodeRight,2,
	       oi_recv->iphif_arr,oi_send->nr,MPI_INT,NodeLeft ,2,
	       MPI_COMM_WORLD,stat);
  MPI_Sendrecv(oi_send->icth0_arr,oi_send->nr,MPI_INT,NodeRight,3,
	       oi_recv->icth0_arr,oi_send->nr,MPI_INT,NodeLeft ,3,
	       MPI_COMM_WORLD,stat);
  MPI_Sendrecv(oi_send->icthf_arr,oi_send->nr,MPI_INT,NodeRight,4,
	       oi_recv->icthf_arr,oi_send->nr,MPI_INT,NodeLeft ,4,
	       MPI_COMM_WORLD,stat);
  MPI_Sendrecv(oi_send->num_pix,oi_send->nr,MPI_INT,NodeRight,5,
	       oi_recv->num_pix,oi_send->nr,MPI_INT,NodeLeft ,5,
	       MPI_COMM_WORLD,stat);
#endif //_HAVE_MPI
}

static void slices2beams(ParamCoLoRe *par,int task)
{
  flouble ***beams=NULL;
  flouble **slices=NULL;
  int ir;

#ifdef _DEBUG
  print_info("Slices to beams\n");
  if(NodeThis==0) timer(0);
#endif //_DEBUG

  if(task==CPY_TASK_DENS) {
    beams=par->dens_beams;
    slices=par->dens_slices;
  }
  else if(task==CPY_TASK_VRAD) {
    beams=par->vrad_beams;
    slices=par->vrad_slices;
  }
  else if(task==CPY_TASK_P_XX) {
    beams=par->p_xx_beams;
    slices=par->p_xx_slices;
  }
  else if(task==CPY_TASK_P_XY) {
    beams=par->p_xy_beams;
    slices=par->p_xy_slices;
  }
  else if(task==CPY_TASK_P_YY) {
    beams=par->p_yy_beams;
    slices=par->p_yy_slices;
  }
  else
    report_error(1,"Wrong task %d\n",task);

  int inode;
#ifdef _HAVE_MPI
  MPI_Status stat;
#endif //_HAVE_MPI

  for(inode=0;inode<NNodes;inode++) {
    //Receive onion_slices from left into onion_sl_dum
    communicate_onion_slices(par->oi_slices,par->oi_sl_dum,&stat);
    for(ir=0;ir<par->oi_slices->nr;ir++) {
      int ib;
      //cth range
      int icth0_sl=par->oi_sl_dum->icth0_arr[ir];
      int icthf_sl=par->oi_sl_dum->icthf_arr[ir];
      int nphi_sl=2*par->oi_slices->nside_arr[ir];

      //Receive slices from left into dum_slices
      flouble *dum_slice=my_malloc(par->oi_sl_dum->num_pix[ir]*sizeof(flouble));
#ifdef _HAVE_MPI
      MPI_Sendrecv(slices[ir],par->oi_slices->num_pix[ir],FLOUBLE_MPI,NodeRight,ir,
      		   dum_slice,par->oi_sl_dum->num_pix[ir],FLOUBLE_MPI,NodeLeft,ir,
		   MPI_COMM_WORLD,&stat);
#else //_HAVE_MPI
      memcpy(dum_slice,slices[ir],par->oi_slices->num_pix[ir]*sizeof(flouble));
#endif //_HAVE_MPI

      //Add dum_slice into beams
      if(par->oi_sl_dum->num_pix[ir]>0) {
	for(ib=0;ib<par->n_beams_here;ib++) {
	  int icth;
	  int icth0_bm=par->oi_beams[ib]->icth0_arr[ir];
	  int icthf_bm=par->oi_beams[ib]->icthf_arr[ir];
	  int iphi0_bm=par->oi_beams[ib]->iphi0_arr[ir];
	  int nphi_bm=par->oi_beams[ib]->iphif_arr[ir]-par->oi_beams[ib]->iphi0_arr[ir]+1;
	  int icth0=MAX(icth0_bm,icth0_sl);
	  int icthf=MIN(icthf_bm,icthf_sl);
	  for(icth=icth0;icth<=icthf;icth++) {
	    int iphi;
	    int ind_cth_bm=(icth-icth0_bm)*nphi_bm;
	    int ind_cth_sl=(icth-icth0_sl)*nphi_sl;
	    for(iphi=0;iphi<nphi_bm;iphi++) {
	      beams[ib][ir][ind_cth_bm+iphi]+=dum_slice[ind_cth_sl+iphi0_bm+iphi];
	    }
	  }
	}
      }

      //Copy dum_slices into current slices
      free(slices[ir]);
      slices[ir]=my_malloc(par->oi_sl_dum->num_pix[ir]*sizeof(flouble));
      memcpy(slices[ir],dum_slice,par->oi_sl_dum->num_pix[ir]*sizeof(flouble));
      free(dum_slice);
    }
    //Copy onion_sl_dum into onion_slices
    memcpy(par->oi_slices->iphi0_arr,par->oi_sl_dum->iphi0_arr,par->oi_slices->nr*sizeof(int));
    memcpy(par->oi_slices->iphif_arr,par->oi_sl_dum->iphif_arr,par->oi_slices->nr*sizeof(int));
    memcpy(par->oi_slices->icth0_arr,par->oi_sl_dum->icth0_arr,par->oi_slices->nr*sizeof(int));
    memcpy(par->oi_slices->icthf_arr,par->oi_sl_dum->icthf_arr,par->oi_slices->nr*sizeof(int));
    memcpy(par->oi_slices->num_pix,par->oi_sl_dum->num_pix,par->oi_slices->nr*sizeof(int));
  }

#ifdef _DEBUG
  if(NodeThis==0) timer(2);
#endif //_DEBUG
}

void pixelize(ParamCoLoRe *par)
{
  //Interpolate into slices
  compute_slices(par,CPY_TASK_DENS);
  compute_slices(par,CPY_TASK_VRAD);
  if(par->do_lensing) {
    compute_slices(par,CPY_TASK_P_XX);
    compute_slices(par,CPY_TASK_P_XY);
    compute_slices(par,CPY_TASK_P_YY);
  }

  end_fftw(par);

  //Communicate into beams
  alloc_beams(par);
  slices2beams(par,CPY_TASK_DENS);
  slices2beams(par,CPY_TASK_VRAD);
  if(par->do_lensing) {
    slices2beams(par,CPY_TASK_P_XX);
    slices2beams(par,CPY_TASK_P_XY);
    slices2beams(par,CPY_TASK_P_YY);
  }

  free_slices(par);
}
