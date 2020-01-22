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

void shear_set_cartesian(ParamCoLoRe *par)
{
  return;
}

void shear_get_local_properties(ParamCoLoRe *par)
{
  return;
}

void shear_distribute(ParamCoLoRe *par)
{
  return;
}

void shear_beams_preproc(ParamCoLoRe *par)
{
  //Sort radii in ascending order
  int ir;
  flouble *r0_i=my_malloc(par->smap->nr*sizeof(flouble));
  flouble *rf_i=my_malloc(par->smap->nr*sizeof(flouble));

  memcpy(r0_i,par->smap->r0,par->smap->nr*sizeof(flouble));
  memcpy(rf_i,par->smap->rf,par->smap->nr*sizeof(flouble));

  int *i_sorted=ind_sort(par->smap->nr,par->smap->r0);
  for(ir=0;ir<par->smap->nr;ir++) {
    par->smap->r0[ir]=r0_i[i_sorted[ir]];
    par->smap->rf[ir]=rf_i[i_sorted[ir]];
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
    for(ipp=0;ipp<2*par->smap->num_pix*par->smap->nr;ipp++) {
      par->smap->data[ipp]=0;
    } //end omp for

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ipp=0;ipp<par->smap->num_pix*par->smap->nr;ipp++) {
      par->smap->nadd[ipp]=1;
    } //end omp for
  } //end omp parallel

  return;
}

void shear_get_beam_properties(ParamCoLoRe *par)
{
  HealpixShells *smap=par->smap;

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,smap)
#endif //_HAVE_OMP
  {
    double idx=par->n_grid/par->l_box;

    //Setup radial decomposition
    int i_r,nr;
    double idr,dr;
    get_radial_params(par->r_max,par->n_grid,&nr,&dr);
    idr=1./dr;

    //Compute index of each source plane
    int *i_r_max_arr=my_malloc(smap->nr*sizeof(int));
    int *i_r_min_arr=my_malloc(smap->nr*sizeof(int));
    double *inv_r_max=my_malloc(smap->nr*sizeof(double));
    for(i_r=0;i_r<smap->nr;i_r++) {
      int i_r_here=(int)(smap->rf[i_r]*idr+0.5);
      inv_r_max[i_r]=1./(i_r_here*dr);
      i_r_max_arr[i_r]=MIN(i_r_here,nr-1);
    }

    i_r_min_arr[0]=0;
    for(i_r=1;i_r<smap->nr;i_r++)
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
    for(ip=0;ip<smap->num_pix;ip++) {
      int ax,added;
      flouble t[6];
      double r1[6],r2[6],xn[3];
      double shear1_1=0,shear1_2=0;
      double shear2_1=0,shear2_2=0;
      double *u=&(smap->pos[3*ip]);
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

      r1[0]=(cth_h*cth_h*cph_h*cph_h-sph_h*sph_h)*prefac;
      r1[1]=(2*cph_h*sph_h*(cth_h*cth_h+1))*prefac;
      r1[2]=(-2*cth_h*sth_h*cph_h)*prefac;
      r1[3]=(cth_h*cth_h*sph_h*sph_h-cph_h*cph_h)*prefac;
      r1[4]=(-2*cth_h*sth_h*sph_h)*prefac;
      r1[5]=(sth_h*sth_h)*prefac;

      r2[0]=(-2*cth_h*cph_h*sph_h)*prefac;
      r2[1]=(2*cth_h*(cph_h*cph_h-sph_h*sph_h))*prefac;
      r2[2]=(2*sth_h*sph_h)*prefac;
      r2[3]=(2*cth_h*sph_h*cph_h)*prefac;
      r2[4]=(-2*sth_h*cph_h)*prefac;
      r2[5]=0;
      for(i_r=0;i_r<smap->nr;i_r++) {
	int irr;
	int irmin=i_r_min_arr[i_r];
	int irmax=i_r_max_arr[i_r];
	for(irr=irmin;irr<=irmax;irr++) {
	  double rm=(irr+0.5)*dr;
	  for(ax=0;ax<3;ax++)
	    xn[ax]=(rm*u[ax]+par->pos_obs[ax])*idx;
	  added=interpolate_from_grid(par,xn,NULL,NULL,t,NULL,NULL,RETURN_TID,INTERP_TYPE_SHEAR);
	  if(added) {
	    double dotp1=0,dotp2=0;
	    for(ax=0;ax<6;ax++) {
	      dotp1+=r1[ax]*t[ax];
	      dotp2+=r2[ax]*t[ax];
            }
	    shear1_1+=dotp1*fac_r_1[irr];
	    shear1_2+=dotp1*fac_r_2[irr];
	    shear2_1+=dotp2*fac_r_1[irr];
	    shear2_2+=dotp2*fac_r_2[irr];
	  }
	}
	smap->data[2*(i_r*smap->num_pix+ip)+0]+=(shear1_1-inv_r_max[i_r]*shear1_2);
	smap->data[2*(i_r*smap->num_pix+ip)+1]+=(shear2_1-inv_r_max[i_r]*shear2_2);
      }
    } //end omp for

#ifdef _HAVE_OMP
#pragma omp single
#endif //_HAVE_OMP
    {
      for(i_r=0;i_r<smap->nr;i_r++) {
        smap->r0[i_r]=1./inv_r_max[i_r];
        smap->rf[i_r]=1./inv_r_max[i_r];
      }
    }

    free(fac_r_1);
    free(fac_r_2);
    free(i_r_max_arr);
    free(i_r_min_arr);
    free(inv_r_max);
  } //end omp parallel

  return;
}

void shear_beams_postproc(ParamCoLoRe *par)
{
  //Set rf to end of cell
  return;
}
