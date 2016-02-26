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
#ifdef _HAVE_FITS
#include <fitsio.h>
#endif //_HAVE_FITS
#ifdef _HAVE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif //_HAVE_HDF5

static ParamCoLoRe *param_colore_new(void)
{
  ParamCoLoRe *par=my_malloc(sizeof(ParamCoLoRe));

  sprintf(par->fnamePk,"default");
  sprintf(par->fnameBz,"default");
  sprintf(par->fnameNz,"default");
  par->output_density=0;
  par->OmegaM=0.3;
  par->OmegaL=0.7;
  par->OmegaB=0.05;
  par->hhub=0.7;
  par->weos=-1.;
  par->n_scal=0.96;
  par->sig8=0.83;
  par->fgrowth_0=-1;
  par->hubble_0=-1;
  par->z_max=1.5;
  par->z_min=0.5;
  par->r_max=-1;
  par->r_min=-1;
  par->r2_smooth=2.0;
  par->do_smoothing=1;
  par->numk=0;
  par->logkmax=1;
  par->logkmin=-3;
  par->idlogk=100;
  par->logkarr=NULL;
  par->pkarr=NULL;
  par->glob_idr=-1;
  par->seed_rng=1234;
  par->z_max=1.;
  par->z_min=0.1;
  par->n_grid=512;
  par->l_box=-1;
  par->nz_here=512;
  par->iz0_here=0;
  sprintf(par->prefixOut,"default");
  par->output_format=0;
  par->grid_dens_f=NULL;
  par->grid_dens=NULL;
  par->grid_vpot_f=NULL;
  par->grid_vpot=NULL;
  par->grid_rvel=NULL;
  par->sigma2_gauss=-1;
  par->nsources=NULL;
  par->spline_bz=NULL;
  par->intacc_bz=NULL;
  par->spline_nz=NULL;
  par->intacc_nz=NULL;
  par->gals=NULL;

  return par;
}

ParamCoLoRe *read_run_params(char *fname)
{
  //////
  // Reads and checks the parameter file
  FILE *fi;
  int n_lin,ii;
  ParamCoLoRe *par=param_colore_new();
  
  //Read parameters from file
  print_info("*** Reading run parameters \n");
  fi=fopen(fname,"r");
  if(fi==NULL) error_open_file(fname);
  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      error_read_line(fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      error_read_line(fname,ii+1);

    if(!strcmp(s1,"prefix_out="))
      sprintf(par->prefixOut,"%s",s2);
    else if(!strcmp(s1,"output_format=")) {
      if(!strcmp(s2,"HDF5")) {
#ifdef _HAVE_HDF5
	par->output_format=2;
#else //_HAVE_HDF5
	report_error(1,"HDF5 format not supported\n");
#endif //_HAVE_HDF5
      }
      else if(!strcmp(s2,"FITS")) {
#ifdef _HAVE_FITS
	par->output_format=1;
#else //_HAVE_FITS
	report_error(1,"FITS format not supported\n");
#endif //_HAVE_FITS
      }
      else if(!strcmp(s2,"ASCII"))
	par->output_format=0;
      else
	report_error(1,"Unrecognized format %s\n",s2);
    }
    else if(!strcmp(s1,"nz_filename="))
      sprintf(par->fnameNz,"%s",s2);
    else if(!strcmp(s1,"bias_filename="))
      sprintf(par->fnameBz,"%s",s2);
    else if(!strcmp(s1,"pk_filename="))
      sprintf(par->fnamePk,"%s",s2);
    else if(!strcmp(s1,"omega_M="))
      par->OmegaM=atof(s2);
    else if(!strcmp(s1,"omega_L="))
      par->OmegaL=atof(s2);
    else if(!strcmp(s1,"omega_B="))
      par->OmegaB=atof(s2);
    else if(!strcmp(s1,"h="))
      par->hhub=atof(s2);
    else if(!strcmp(s1,"w="))
      par->weos=atof(s2);
    else if(!strcmp(s1,"ns="))
      par->n_scal=atof(s2);
    else if(!strcmp(s1,"sigma_8="))
      par->sig8=atof(s2);
    else if(!strcmp(s1,"r_smooth="))
      par->r2_smooth=atof(s2);
    else if(!strcmp(s1,"z_min="))
      par->z_min=atof(s2);
    else if(!strcmp(s1,"z_max="))
      par->z_max=atof(s2);
    else if(!strcmp(s1,"n_grid="))
      par->n_grid=atoi(s2);
    else if(!strcmp(s1,"seed="))
      par->seed_rng=atoi(s2);
    else if(!strcmp(s1,"output_density="))
      par->output_density=atoi(s2);
    else
      fprintf(stderr,"CoLoRe: Unknown parameter %s\n",s1);
  }
  fclose(fi);

  if(par->r2_smooth>0) {
    par->r2_smooth=pow(par->r2_smooth,2);
    par->do_smoothing=1;
  }
  else
    par->do_smoothing=0;
  cosmo_set(par);
  init_fftw(par);

  double dk=2*M_PI/par->l_box;
  print_info("Run parameters: \n");
  print_info("  %.3lf < z < %.3lf\n",par->z_min,par->z_max);
  print_info("  %.3lf < r/(Mpc/h) < %.3lf\n",par->r_min,par->r_max);
  print_info("  L_box = %.3lf Mpc/h, N_grid = %d \n",par->l_box,par->n_grid);
  print_info("  Scales resolved: %.3lE < k < %.3lE h/Mpc\n",dk,0.5*(par->n_grid-1)*dk);
  print_info("  Fourier-space resolution: dk = %.3lE h/Mpc\n",dk);
  print_info("  Real-space resolution: dx = %.3lE Mpc/h\n",par->l_box/par->n_grid);
  if(par->do_smoothing)
    print_info("  Density field pre-smoothed on scales: x_s = %.3lE Mpc/h\n",sqrt(par->r2_smooth));
  else
    print_info("  No extra smoothing\n");
  print_info("\n");
  
  return par;
}

void write_grid(ParamCoLoRe *par)
{
  FILE *fo;
  char fname[256];
  int iz;
  int ngx=2*(par->n_grid/2+1);
  int size_flouble=sizeof(flouble);

  if(NodeThis==0) timer(0);
  print_info("*** Writing density field (native format)\n");
  sprintf(fname,"%s_dens_%d.dat",par->prefixOut,NodeThis);
  fo=fopen(fname,"wb");
  if(fo==NULL) error_open_file(fname);
  my_fwrite(&NNodes,sizeof(int),1,fo);
  my_fwrite(&size_flouble,sizeof(int),1,fo);
  my_fwrite(&(par->l_box),sizeof(double),1,fo);
  my_fwrite(&(par->n_grid),sizeof(int),1,fo);
  my_fwrite(&(par->nz_here),sizeof(int),1,fo);
  my_fwrite(&(par->iz0_here),sizeof(int),1,fo);
  for(iz=0;iz<par->nz_here;iz++) {
    int iy;
    for(iy=0;iy<par->n_grid;iy++) {
      lint index0=ngx*((lint)(iy+iz*par->n_grid));
      my_fwrite(&(par->grid_dens[index0]),sizeof(flouble),par->n_grid,fo);
    }
  }
  fclose(fo);
  if(NodeThis==0) timer(2);
  print_info("\n");
}

void write_catalog(ParamCoLoRe *par)
{
  char fname[256];
  
  if(NodeThis==0) timer(0);
  if(par->output_format==2) { //HDF5
#ifdef _HAVE_HDF5
    hid_t file_id,gal_types[4];
    hsize_t chunk_size=128;
    size_t dst_offset[4]={HOFFSET(Gal,ra),HOFFSET(Gal,dec),HOFFSET(Gal,z0),HOFFSET(Gal,dz_rsd)};
    const char *names[4]={"RA" ,"DEC","Z_COSMO","DZ_RSD"};
    char *tunit[4]=      {"DEG","DEG","NA"     ,"NA"    };
    gal_types[0]=H5T_NATIVE_FLOAT;
    gal_types[1]=H5T_NATIVE_FLOAT;
    gal_types[2]=H5T_NATIVE_FLOAT;
    gal_types[3]=H5T_NATIVE_FLOAT;

    print_info("*** Writing output (HDF5)\n");
    sprintf(fname,"%s_%d.h5",par->prefixOut,NodeThis);

    //Create file
    file_id=H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    //Write table
    H5TBmake_table("source_data",file_id,"/sources",4,par->nsources_this,sizeof(Gal),
		   names,dst_offset,gal_types,chunk_size,NULL,0,par->gals);
    H5LTset_attribute_string(file_id,"/sources","FIELD_0_UNITS",tunit[0]);
    H5LTset_attribute_string(file_id,"/sources","FIELD_1_UNITS",tunit[1]);
    H5LTset_attribute_string(file_id,"/sources","FIELD_2_UNITS",tunit[2]);
    H5LTset_attribute_string(file_id,"/sources","FIELD_3_UNITS",tunit[3]);

    //End file
    H5Fclose(file_id);
#else //_HAVE_HDF5
    printf("HDF5 not supported\n");
#endif //_HAVE_HDF5
  }
  else if(par->output_format==1) { //FITS
#ifdef _HAVE_FITS
    long ii,row_here=0,nrw=0;
    int status=0;
    fitsfile *fptr;
    float *ra_arr,*dec_arr,*z0_arr,*rsd_arr;
    char *ttype[]={"RA" ,"DEC","Z_COSMO","DZ_RSD"};
    char *tform[]={"1E" ,"1E" ,"1E"     ,"1E"    };
    char *tunit[]={"DEG","DEG","NA"     ,"NA"    };

    print_info("*** Writing output (FITS)\n");
    sprintf(fname,"!%s_%d.fits",par->prefixOut,NodeThis);

    fits_create_file(&fptr,fname,&status);
    fits_create_tbl(fptr,BINARY_TBL,0,4,ttype,tform,tunit,NULL,&status);

    fits_get_rowsize(fptr,&nrw,&status);
    ra_arr=my_malloc(nrw*sizeof(float));
    dec_arr=my_malloc(nrw*sizeof(float));
    z0_arr=my_malloc(nrw*sizeof(float));
    rsd_arr=my_malloc(nrw*sizeof(float));
    while(row_here<par->nsources_this) {
      long nrw_here;
      if(row_here+nrw>par->nsources_this)
	nrw_here=par->nsources_this-row_here;
      else
	nrw_here=nrw;

      for(ii=0;ii<nrw_here;ii++) {
	ra_arr[ii]=par->gals[row_here+ii].ra;
	dec_arr[ii]=par->gals[row_here+ii].dec;
	z0_arr[ii]=par->gals[row_here+ii].z0;
	rsd_arr[ii]=par->gals[row_here+ii].dz_rsd;
      }
      fits_write_col(fptr,TFLOAT,1,row_here+1,1,nrw_here,ra_arr,&status);
      fits_write_col(fptr,TFLOAT,2,row_here+1,1,nrw_here,dec_arr,&status);
      fits_write_col(fptr,TFLOAT,3,row_here+1,1,nrw_here,z0_arr,&status);
      fits_write_col(fptr,TFLOAT,4,row_here+1,1,nrw_here,rsd_arr,&status);
    
      row_here+=nrw_here;
    }
    fits_close_file(fptr,&status);
    free(ra_arr);
    free(dec_arr);
    free(z0_arr);
    free(rsd_arr);
#else //_HAVE_FITS
    printf("FITS not supported\n");
#endif //_HAVE_FITS
  }
  else   {
    print_info("*** Writing output (ASCII)\n");
    sprintf(fname,"%s_%d.txt",par->prefixOut,NodeThis);

    lint jj;
    FILE *fil=fopen(fname,"w");
    if(fil==NULL) error_open_file(fname);
    fprintf(fil,"#[1]RA, [2]dec, [3]z0, [4]dz_RSD\n");
    for(jj=0;jj<par->nsources_this;jj++) {
      fprintf(fil,"%E %E %E %E\n",
	      par->gals[jj].ra,par->gals[jj].dec,
	      par->gals[jj].z0,par->gals[jj].dz_rsd);
    }
    fclose(fil);
  }
  if(NodeThis==0) timer(2);
  print_info("\n");
}

void param_colore_free(ParamCoLoRe *par)
{
  free(par->logkarr);
  free(par->pkarr);
#ifdef _SPREC
  if(par->grid_dens_f!=NULL)
    fftwf_free(par->grid_dens_f);
  fftwf_free(par->grid_vpot_f);
#else //_SPREC
  if(par->grid_dens_f!=NULL)
    fftw_free(par->grid_dens_f);
  fftw_free(par->grid_vpot_f);
#endif //_SPREC
#ifdef _HAVE_MPI
  free(par->slice_left);
  free(par->slice_right);
#endif //_HAVE_MPI
  free(par->grid_rvel);
  free(par->nsources);
  gsl_spline_free(par->spline_bz);
  gsl_spline_free(par->spline_nz);
  gsl_interp_accel_free(par->intacc_bz);
  gsl_interp_accel_free(par->intacc_nz);
  if(par->gals!=NULL)
    free(par->gals);
  end_fftw();
}
