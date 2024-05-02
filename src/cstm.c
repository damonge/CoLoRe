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

void cstm_set_cartesian(ParamCoLoRe *par)
{
  return;
}

void cstm_get_local_properties(ParamCoLoRe *par)
{
  return;
}

void cstm_distribute(ParamCoLoRe *par)
{
  return;
}

static void cstm_beams_preproc_single(ParamCoLoRe *par, int ipop)
{
  //Zero all data and clear
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,ipop)
#endif //_HAVE_OMP
  {
    long ipp;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ipp=0;ipp<par->cstm[ipop]->num_pix;ipp++) {
      par->cstm[ipop]->data[ipp]=0;
      par->cstm[ipop]->nadd[ipp]=1;
    } //end omp for
  } //end omp parallel

  return;
}

void cstm_beams_preproc(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_cstm;ipop++)
    cstm_beams_preproc_single(par,ipop);
}

static void cstm_get_beam_properties_single(ParamCoLoRe *par,int ipop)
{
  HealpixShells *cmap=par->cstm[ipop];

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,cmap,ipop)
#endif //_HAVE_OMP
  {
    double idx=par->n_grid/par->l_box;

    //Setup radial decomposition
    int ir,nr;
    double dr;
    get_radial_params(par->r_max,par->n_grid,&nr,&dr);

    //Precompute kernel, bias and normalization
    double k_max=-1E100;
    double *kz=my_malloc(nr*sizeof(double));
    double *bz=my_malloc(nr*sizeof(double));
    double *normz=my_malloc(nr*sizeof(double));
    for(ir=0;ir<nr;ir++) {
      double rm=(ir+0.5)*dr;
      //Kernel * H(z)
      kz[ir]=get_bg(par,rm,BG_KZ_CSTM,ipop)/get_bg(par,rm,BG_IH,0);
      bz[ir]=get_bg(par,rm,BG_BZ_CSTM,ipop);
      normz[ir]=get_bg(par,rm,BG_NORM_CSTM,ipop);
      if(fabs(kz[ir])>k_max)
	k_max=fabs(kz[ir]);
    }
    //Find integration edges
    int ir_min=0, ir_max=nr-1;

    for(ir=0;ir<nr;ir++) {
      if(fabs(kz[ir])>1E-4*k_max) {
	ir_min=ir;
	break;
      }
    }
    for(ir=nr-1;ir>=0;ir--) {
      if(fabs(kz[ir])>1E-4*k_max) {
	ir_max=ir;
	break;
      }
    }
    if(ir_max < ir_min)
      report_error(1, "Custom kernel has no suppport");


    long ip;
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ip=0;ip<cmap->num_pix;ip++) {
      int ax,irr,added;
      flouble d;
      double xn[3];
      double *u=&(cmap->pos[3*ip]);
      double cval=0;

      for(irr=ir_min;irr<=ir_max;irr++) {
	double rm=(irr+0.5)*dr;
	for(ax=0;ax<3;ax++)
	  xn[ax]=(rm*u[ax]+par->pos_obs[ax])*idx;
	added=interpolate_from_grid(par,xn,&d,NULL,NULL,NULL,NULL,RETURN_DENS,INTERP_CIC);
	if(added)
	  cval+=kz[irr]*(bias_model(d,bz[irr], -1)*normz[irr]-1);
      }
      cmap->data[ip]+=cval*dr;
    } //end omp for

    free(kz);
    free(bz);
    free(normz);
  } //end omp parallel

  return;
}

void cstm_get_beam_properties(ParamCoLoRe *par)
{
  int ipop;
  for(ipop=0;ipop<par->n_cstm;ipop++)
    cstm_get_beam_properties_single(par,ipop);
}

void cstm_beams_postproc(ParamCoLoRe *par)
{
  //Set rf to end of cell
  return;
}
