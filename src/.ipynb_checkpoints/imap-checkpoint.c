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

static int get_r_index_imap(HealpixShells *sh,double r,int ir_start)
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
    if((ir0==-1) || (ir0==sh->nr))
      gotit=1;
    else {
      if(r<sh->r0[ir0])
	ir0--;
      else {
	if(r>=sh->rf[ir0])
	  ir0++;
	else
	  gotit=1;
      }
    }
  }
  
  return ir0;
}

static void imap_preproc(HealpixShells *imap)
{
  int ir;
  flouble *r0_i=my_malloc(imap->nr*sizeof(flouble));
  flouble *rf_i=my_malloc(imap->nr*sizeof(flouble));

  memcpy(r0_i,imap->r0,imap->nr*sizeof(flouble));
  memcpy(rf_i,imap->rf,imap->nr*sizeof(flouble));

  int *i_sorted=ind_sort(imap->nr,imap->r0);
  for(ir=0;ir<imap->nr;ir++) {
    imap->r0[ir]=r0_i[i_sorted[ir]];
    imap->rf[ir]=rf_i[i_sorted[ir]];
  }
  free(r0_i);
  free(rf_i);
  free(i_sorted);
  
  //All nodes must store the entire sky
  imap->num_pix=he_nside2npix(imap->nside);
  free(imap->listpix);
  imap->listpix=my_malloc(sizeof(long)); //No need to initialize this really
  free(imap->pos);
  imap->pos=my_malloc(sizeof(double)); //No need to initialize thi either
  free(imap->data);
  imap->data=my_calloc(imap->nr*imap->num_pix,sizeof(flouble));
  free(imap->nadd);
  imap->nadd=my_calloc(imap->nr*imap->num_pix,sizeof(flouble));
}

static void imap_set_cartesian_single(ParamCoLoRe *par,int ipop)
{
  if(NodeThis==0) timer(0);
  HealpixShells *imap=par->imap[ipop];
  imap_preproc(imap);

  print_info(" %d-th intensity mapping species\n",ipop);
  if(NodeThis==0) timer(0);
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,imap,ipop)
#endif //_HAVE_OMP
  {
    long iz;
    double dx=par->l_box/par->n_grid;
    int ngx=2*(par->n_grid/2+1);
    double factor_vel=-par->fgrowth_0/(1.5*par->hubble_0*par->OmegaM);

    double pixel_size=sqrt(4*M_PI/he_nside2npix(imap->nside));
    double rmin_here=imap->r0[0]-20.;
    double rmax_here=imap->rf[imap->nr-1]+20.;
    int nsub_lo,nsub_hi;
    int *nsub_array=my_malloc(imap->nr*sizeof(int));
    for(iz=0;iz<imap->nr;iz++) {
      double r0=imap->r0[iz];
      double rf=imap->rf[iz];
      double vol=(rf*rf*rf-r0*r0*r0)*pixel_size*pixel_size/3;
      double sct=r0*pixel_size; //Scale associated to pixel size
      double scr=rf-r0; //Scale associated to the width of the shell
      double scv=pow(vol,0.333333); //Scale associated to the overall volume
      double sc0=fmin(scv,fmin(sct,scr)); //We now pick the minimum of the previous 3

      nsub_array[iz]=(int)(dx/sc0+0.5)+1;
      //#ifdef _HAVE_OMP
      //      if(omp_get_thread_num()==0)
      //#endif //_HAVE_OMP
      //	print_info("%ld %lf %lf %d %lf %lf\n",iz,r0,rf,nsub_array[iz],rmin_here,rmax_here);
    }
    nsub_lo=nsub_array[0];
    nsub_hi=nsub_array[imap->nr-1];

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      int irad=0;
      long indexz=iz*((long)(ngx*par->n_grid));
      double z0=(iz+par->iz0_here)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	double y0=iy*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  double x0=ix*dx-par->pos_obs[0];
	  long index=ix+indexy+indexz;
	  double r0=sqrt(x0*x0+y0*y0+z0*z0);
	  double tmean=get_bg(par,r0,BG_TZ_IMAP,ipop);
	  if((r0<=rmax_here) && (r0>=rmin_here)) {
	    if(tmean>0) {
	      int nsub_here,izz;
	      double dx_sub;
	      double bias=get_bg(par,r0,BG_BZ_IMAP,ipop);
	      double dnorm=get_bg(par,r0,BG_NORM_IMAP,ipop);
	      double rvel=factor_vel*get_rvel(par,ix,iy,iz,x0,y0,z0,r0);
	      double dr_rsd=rvel*get_bg(par,r0,BG_V1,0)*get_bg(par,r0,BG_IH,0);
	      double temp=tmean*bias_model(par->grid_dens[index],bias)*dnorm;
	      irad=get_r_index_imap(imap,r0,irad);
	      if(irad<0)
		nsub_here=nsub_lo;
	      else if(irad>=imap->nr)
		nsub_here=nsub_hi;
	      else
		nsub_here=nsub_array[irad];
	      dx_sub=dx/nsub_here;
	      for(izz=0;izz<nsub_here;izz++) {
		int iyy;
		double z=z0+(izz+0.5)*dx_sub;
		for(iyy=0;iyy<nsub_here;iyy++) {
		  int ixx;
		  double y=y0+(iyy+0.5)*dx_sub;
		  for(ixx=0;ixx<nsub_here;ixx++) {
		    double cth,phi,r;
		    double x=x0+(ixx+0.5)*dx_sub;
		    cart2sph(x,y,z,&r,&cth,&phi);
		    irad=get_r_index_imap(imap,r+dr_rsd,irad);
		    if((irad>=0) && (irad<imap->nr)) {
		      long irad_t=irad*imap->num_pix;
		      long pix_id=he_ang2pix(imap->nside,cth,phi);
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
		      imap->data[irad_t+pix_id]+=temp;
#ifdef _HAVE_OMP
#pragma omp atomic
#endif //_HAVE_OMP
		      imap->nadd[irad_t+pix_id]++;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } //end omp for

    free(nsub_array);
  } //end omp parallel
  if(NodeThis==0) timer(2);
}

void imap_set_cartesian(ParamCoLoRe *par)
{
  int ipop;

  //First, compute lensing Hessian
  print_info("*** Filling up intensity maps\n");
  for(ipop=0;ipop<par->n_imap;ipop++)
    imap_set_cartesian_single(par,ipop);
  print_info("\n");
}

static void imap_distribute_single(ParamCoLoRe *par,int ipop)
{
  return;
}

void imap_distribute(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_imap;ipop++)
    imap_distribute_single(par,ipop);
}

static void imap_get_local_properties_single(ParamCoLoRe *par,int ipop)
{
  return;
}

void imap_get_local_properties(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_imap;ipop++)
    imap_get_local_properties_single(par,ipop);
}

static void imap_beams_preproc_single(ParamCoLoRe *par,int ipop)
{
  return;
}

void imap_beams_preproc(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_imap;ipop++)
    imap_beams_preproc_single(par,ipop);
}

static void imap_get_beam_properties_single(ParamCoLoRe *par,int ipop)
{
  return;
}

void imap_get_beam_properties(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_imap;ipop++)
    imap_get_beam_properties_single(par,ipop);
}

static void imap_beams_postproc_single(ParamCoLoRe *par,int ipop)
{
  return;
}

void imap_beams_postproc(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_imap;ipop++)
    imap_beams_postproc_single(par,ipop);
}
