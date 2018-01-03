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

void kappa_set_cartesian(ParamCoLoRe *par)
{
  return;
}

void kappa_get_local_properties(ParamCoLoRe *par)
{
  return;
}

void kappa_distribute(ParamCoLoRe *par)
{
  return;
}

void kappa_beams_preproc(ParamCoLoRe *par)
{
  //Sort radii in ascending order
  int ir;
  flouble *r0_i=my_malloc(par->kmap->nr*sizeof(flouble));
  flouble *rf_i=my_malloc(par->kmap->nr*sizeof(flouble));

  memcpy(r0_i,par->kmap->r0,par->kmap->nr*sizeof(flouble));
  memcpy(rf_i,par->kmap->rf,par->kmap->nr*sizeof(flouble));

  int *i_sorted=ind_sort(par->kmap->nr,par->kmap->r0);
  for(ir=0;ir<par->kmap->nr;ir++) {
    par->kmap->r0[ir]=r0_i[i_sorted[ir]];
    par->kmap->rf[ir]=rf_i[i_sorted[ir]];
  }
  free(r0_i);
  free(rf_i);
  free(i_sorted);

  //Zero all data and clear
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par)
#endif //_HAVE_OMP
  {
    long ipp;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ipp=0;ipp<par->kmap->num_pix*par->kmap->nr;ipp++) {
      par->kmap->data[ipp]=0;
      par->kmap->nadd[ipp]=1;
    } //end omp for
  } //end omp parallel

  return;
}

void kappa_get_beam_properties(ParamCoLoRe *par)
{
  HealpixShells *kmap=par->kmap;

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,kmap)
#endif //_HAVE_OMP
  {
    double idx=par->n_grid/par->l_box;

    //Setup radial decomposition
    int i_r,nr;
    double idr,dr;
    get_radial_params(par->r_max,par->n_grid,&nr,&dr);
    idr=1./dr;
  
    //Compute index of each source plane
    int *i_r_max_arr=my_malloc(kmap->nr*sizeof(int));
    int *i_r_min_arr=my_malloc(kmap->nr*sizeof(int));
    double *inv_r_max=my_malloc(kmap->nr*sizeof(double));
    for(i_r=0;i_r<kmap->nr;i_r++) {
      int i_r_here=(int)(kmap->rf[i_r]*idr+0.5);
      inv_r_max[i_r]=1./(i_r_here*dr);
      i_r_max_arr[i_r]=MIN(i_r_here,nr-1);
    }
    i_r_min_arr[0]=0;
    for(i_r=1;i_r<kmap->nr;i_r++)
      i_r_min_arr[i_r]=i_r_max_arr[i_r-1]+1;

    //Fill up integral kernels
    double *fac_r_1=my_malloc(nr*sizeof(double));
    double *fac_r_2=my_malloc(nr*sizeof(double));
    for(i_r=0;i_r<nr;i_r++) {
      double rm=(i_r+0.5)*dr;
      double pg=get_bg(par,rm,BG_D1,0)*(1+get_bg(par,rm,BG_Z,0));
      fac_r_1[i_r]=rm*pg*dr;
      fac_r_2[i_r]=rm*rm*pg*dr;
    }

    long ip;
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ip=0;ip<kmap->num_pix;ip++) {
      int ax,added;
      flouble t[6];
      double rot[6],xn[3];
      double kappa_1=0,kappa_2=0;
      double *u=&(kmap->pos[3*ip]);
      double prefac=idx*idx;
      double cth_h=1,sth_h=0,cph_h=1,sph_h=0;
      
      cth_h=u[2];
      if(cth_h>=1) cth_h=1;
      if(cth_h<=-1) cth_h=-1;
      sth_h=sqrt((1-cth_h)*(1+cth_h));
      if(sth_h!=0) {
	cph_h=u[0]/sth_h;
	sph_h=u[1]/sth_h;
      }
      
      rot[0]=(cth_h*cth_h*cph_h*cph_h+sph_h*sph_h)*prefac;
      rot[1]=(2*cph_h*sph_h*(cth_h*cth_h-1))*prefac;
      rot[2]=(-2*cth_h*sth_h*cph_h)*prefac;
      rot[3]=(cth_h*cth_h*sph_h*sph_h+cph_h*cph_h)*prefac;
      rot[4]=(-2*cth_h*sth_h*sph_h)*prefac;
      rot[5]=(sth_h*sth_h)*prefac;
      for(i_r=0;i_r<kmap->nr;i_r++) {
	int irr;
	int irmin=i_r_min_arr[i_r];
	int irmax=i_r_max_arr[i_r];
	for(irr=irmin;irr<=irmax;irr++) {
	  double rm=(irr+0.5)*dr;
	  for(ax=0;ax<3;ax++)
	    xn[ax]=(rm*u[ax]+par->pos_obs[ax])*idx;
	  added=interpolate_from_grid(par,xn,NULL,NULL,t,NULL,RETURN_TID,INTERP_TYPE_SHEAR);
	  if(added) {
	    double dotp=0;
	    for(ax=0;ax<6;ax++)
	      dotp+=rot[ax]*t[ax];
	    kappa_1+=dotp*fac_r_1[irr];
	    kappa_2+=dotp*fac_r_2[irr];
	  }
	}
	kmap->data[i_r*kmap->num_pix+ip]+=(kappa_1-inv_r_max[i_r]*kappa_2);
      }	
    } //end omp for

    free(fac_r_1);
    free(fac_r_2);
    free(i_r_max_arr);
    free(i_r_min_arr);
    free(inv_r_max);
  } //end omp parallel

  return;
}

void kappa_beams_postproc(ParamCoLoRe *par)
{
  //Set rf to end of cell
  return;
}
