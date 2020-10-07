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

void lensing_set_cartesian(ParamCoLoRe *par)
{
  return;
}

void lensing_get_local_properties(ParamCoLoRe *par)
{
  return;
}

void lensing_distribute(ParamCoLoRe *par)
{
  return;
}

void lensing_beams_preproc(ParamCoLoRe *par)
{
  //Sort radii in ascending order
  int ir;
  flouble *r_i=my_malloc(par->smap->nr*sizeof(flouble));

  memcpy(r_i,par->smap->r,par->smap->nr*sizeof(flouble));

  int *i_sorted=ind_sort(par->smap->nr,par->smap->r);
  for(ir=0;ir<par->smap->nr;ir++) {
    par->smap->r[ir]=r_i[i_sorted[ir]];
  }
  free(r_i);
  free(i_sorted);

  //Zero all data and clear
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par)
#endif //_HAVE_OMP
  {
    long ir, ib, ipp;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ir=0;ir<par->smap->nr;ir++) {
      for(ib=0;ib<par->smap->nbeams;ib++) {
        for(ipp=0;ipp<5*par->smap->num_pix_per_beam[ir];ipp++)
          par->smap->data[ib][ir][ipp]=0;
      }
    } //end omp for
  } //end omp parallel

  return;
}

void lensing_get_beam_properties(ParamCoLoRe *par)
{
  HealpixShellsAdaptive *smap=par->smap;

#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par,smap)
#endif //_HAVE_OMP
  {
    double idx=par->n_grid/par->l_box;
    long npix_hi=smap->num_pix_per_beam[smap->nr-1];

    //Setup radial decomposition
    int i_r,nr;
    double idr,dr;
    get_radial_params(par->r_max,par->n_grid,&nr,&dr);
    idr=1./dr;

    //Compute index of each source plane
    int *i_r_max_arr=my_malloc(smap->nr*sizeof(int));
    int *i_r_min_arr=my_malloc(smap->nr*sizeof(int));
    long *npix_ratio=my_malloc(smap->nr*sizeof(long));
    double *inv_npix_ratio=my_malloc(smap->nr*sizeof(double));
    double *inv_r_max=my_malloc(smap->nr*sizeof(double));
    for(i_r=0;i_r<smap->nr;i_r++) {
      int i_r_here=(int)(smap->r[i_r]*idr+0.5);
      inv_r_max[i_r]=1./(i_r_here*dr);
      i_r_max_arr[i_r]=MIN(i_r_here,nr-1);
      npix_ratio[i_r]=npix_hi/smap->num_pix_per_beam[i_r];
      inv_npix_ratio[i_r]=1./((double)(npix_ratio[i_r]));
    }
  
    i_r_min_arr[0]=0;
    for(i_r=1;i_r<smap->nr;i_r++)
      i_r_min_arr[i_r]=i_r_max_arr[i_r-1]+1;

    //Fill up integral kernels
    double *fac_r_0=my_malloc(nr*sizeof(double));
    double *fac_r_1=my_malloc(nr*sizeof(double));
    double *fac_r_2=my_malloc(nr*sizeof(double));
    for(i_r=0;i_r<nr;i_r++) {
      double rm=(i_r+0.5)*dr;
      double pg=get_bg(par,rm,BG_D1,0)*(1+get_bg(par,rm,BG_Z,0));
      fac_r_0[i_r]=pg*dr;
      fac_r_1[i_r]=rm*pg*dr;
      fac_r_2[i_r]=rm*rm*pg*dr;
    }

    int ib;
    for(ib=0;ib<smap->nbeams;ib++) {
      //We loop over all pixels in the last plane (since this has the
      //highest resolution).
      long ip;
#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
      for(ip=0;ip<npix_hi;ip++) {
        int ax,added;
        flouble t[6],v[3];
        double xn[3], u_x[3], u_y[3];
        double r_k[6], r_e1[6],r_e2[6];
        double dx_0=0,dx_1=0;
        double dy_0=0,dy_1=0;
        double kappa_1=0,kappa_2=0;
        double shear1_1=0,shear1_2=0;
        double shear2_1=0,shear2_2=0;
        double *u=&(smap->pos[ib][3*ip]);
        double prefac=idx*idx; //1/Dx^2
        double prefac_m=0.5*idx; //1/(2 * Dx)
        double cth_h=1,sth_h=0,cph_h=1,sph_h=0;

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
        for(i_r=0;i_r<smap->nr;i_r++) {
          int irr;
          int irmin=i_r_min_arr[i_r];
          int irmax=i_r_max_arr[i_r];
          long ipix_this=ip/npix_ratio[i_r];
          for(irr=irmin;irr<=irmax;irr++) {
            double rm=(irr+0.5)*dr;
            for(ax=0;ax<3;ax++)
              xn[ax]=(rm*u[ax]+par->pos_obs[ax])*idx;
            added=interpolate_from_grid(par,xn,NULL,v,t,NULL,NULL,RETURN_VEL | RETURN_TID,
                                        INTERP_TYPE_LENSING);
            if(added) {
              double dotk=0,dote1=0,dote2=0,dotvx=0,dotvy=0;
              for(ax=0;ax<6;ax++) {
                dote1+=r_e1[ax]*t[ax];
                dote2+=r_e2[ax]*t[ax];
                dotk +=r_k[ax]*t[ax];
              }
              for(ax=0;ax<3;ax++) {
                dotvx+=u_x[ax]*v[ax];
                dotvy+=u_y[ax]*v[ax];
              }
              dx_0+=dotvx*fac_r_0[irr];
              dx_1+=dotvx*fac_r_1[irr];
              dy_0+=dotvy*fac_r_0[irr];
              dy_1+=dotvy*fac_r_1[irr];
              kappa_1+=dotk*fac_r_1[irr];
              kappa_2+=dotk*fac_r_2[irr];
              shear1_1+=dote1*fac_r_1[irr];
              shear1_2+=dote1*fac_r_2[irr];
              shear2_1+=dote2*fac_r_1[irr];
              shear2_2+=dote2*fac_r_2[irr];
            }
          }
          if(ipix_this>=smap->num_pix_per_beam[i_r]) {
            printf("SHIT\n");
            exit(1);
          }
          smap->data[ib][i_r][5*ipix_this+0]+=(shear1_1-inv_r_max[i_r]*shear1_2)*inv_npix_ratio[i_r];
          smap->data[ib][i_r][5*ipix_this+1]+=(shear2_1-inv_r_max[i_r]*shear2_2)*inv_npix_ratio[i_r];
          smap->data[ib][i_r][5*ipix_this+2]+=(kappa_1-inv_r_max[i_r]*kappa_2)*inv_npix_ratio[i_r];
          smap->data[ib][i_r][5*ipix_this+3]+=2*(dx_0-inv_r_max[i_r]*dx_1)*inv_npix_ratio[i_r];
          smap->data[ib][i_r][5*ipix_this+4]+=2*(dy_0-inv_r_max[i_r]*dy_1)*inv_npix_ratio[i_r];
        }
      } //end omp for
    }

#ifdef _HAVE_OMP
#pragma omp single
#endif //_HAVE_OMP
    {
      for(i_r=0;i_r<smap->nr;i_r++)
        smap->r[i_r]=1./inv_r_max[i_r];
    }

    free(fac_r_0);
    free(fac_r_1);
    free(fac_r_2);
    free(i_r_max_arr);
    free(i_r_min_arr);
    free(inv_r_max);
    free(npix_ratio);
    free(inv_npix_ratio);
  } //end omp parallel

  return;
}

void lensing_beams_postproc(ParamCoLoRe *par)
{
  //Set rf to end of cell
  return;
}
