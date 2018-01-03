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

void isw_set_cartesian(ParamCoLoRe *par)
{
  return;
}

void isw_get_local_properties(ParamCoLoRe *par)
{
  return;
}

void isw_distribute(ParamCoLoRe *par)
{
  return;
}

void isw_beams_preproc(ParamCoLoRe *par)
{
  //Sort radii in ascending order
  int ir;
  flouble *r0_i=my_malloc(par->pd_map->nr*sizeof(flouble));
  flouble *rf_i=my_malloc(par->pd_map->nr*sizeof(flouble));

  memcpy(r0_i,par->pd_map->r0,par->pd_map->nr*sizeof(flouble));
  memcpy(rf_i,par->pd_map->rf,par->pd_map->nr*sizeof(flouble));

  int *i_sorted=ind_sort(par->pd_map->nr,par->pd_map->r0);
  for(ir=0;ir<par->pd_map->nr;ir++) {
    par->pd_map->r0[ir]=r0_i[i_sorted[ir]];
    par->pd_map->rf[ir]=rf_i[i_sorted[ir]];
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
    for(ipp=0;ipp<par->pd_map->num_pix*par->pd_map->nr;ipp++) {
      par->pd_map->data[ipp]=0;
      par->pd_map->nadd[ipp]=1;
    } //end omp for
  } //end omp parallel

  return;
}

void isw_get_beam_properties(ParamCoLoRe *par)
{
  HealpixShells *pd_map=par->pd_map;

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,pd_map)
#endif //_HAVE_OMP
  {
    double idx=par->n_grid/par->l_box;

    //Setup radial decomposition
    int i_r,nr;
    double idr,dr;
    get_radial_params(par->r_max,par->n_grid,&nr,&dr);
    idr=1./dr;
  
    //Compute index of each source plane
    int *i_r_max_arr=my_malloc(pd_map->nr*sizeof(int));
    int *i_r_min_arr=my_malloc(pd_map->nr*sizeof(int));
    for(i_r=0;i_r<pd_map->nr;i_r++) {
      int i_r_here=(int)(pd_map->rf[i_r]*idr+0.5);
      i_r_max_arr[i_r]=MIN(i_r_here,nr-1);
    }
    i_r_min_arr[0]=0;
    for(i_r=1;i_r<pd_map->nr;i_r++)
      i_r_min_arr[i_r]=i_r_max_arr[i_r-1]+1;

    //Fill up integral kernels
    double *fac_r=my_malloc(nr*sizeof(double));
    for(i_r=0;i_r<nr;i_r++) {
      double rm=(i_r+0.5)*dr;
      double pdg=get_bg(par,rm,BG_PD,0);
      fac_r[i_r]=2*pdg*dr;
    }

    long ip;
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ip=0;ip<pd_map->num_pix;ip++) {
      int ax,added;
      flouble pd;
      double xn[3];
      double isw=0;
      double *u=&(pd_map->pos[3*ip]);

      for(i_r=0;i_r<pd_map->nr;i_r++) {
	int irr;
	int irmin=i_r_min_arr[i_r];
	int irmax=i_r_max_arr[i_r];
	for(irr=irmin;irr<=irmax;irr++) {
	  double rm=(irr+0.5)*dr;
	  for(ax=0;ax<3;ax++)
	    xn[ax]=(rm*u[ax]+par->pos_obs[ax])*idx;
	  added=interpolate_from_grid(par,xn,NULL,NULL,NULL,&pd,RETURN_PDOT,INTERP_TYPE_SHEAR);
	  if(added)
	    isw+=pd*fac_r[irr];
	}
	pd_map->data[i_r*pd_map->num_pix+ip]+=isw;
      }	
    } //end omp for

    free(fac_r);
    free(i_r_max_arr);
    free(i_r_min_arr);
  } //end omp parallel

  return;
}

void isw_beams_postproc(ParamCoLoRe *par)
{
  //Set rf to end of cell
  return;
}
