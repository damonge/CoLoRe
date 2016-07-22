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

static ParamCoLoRe *param_colore_new(void)
{
  int ii;
  ParamCoLoRe *par=my_malloc(sizeof(ParamCoLoRe));

#ifdef _DEBUG
  par->f_dbg=NULL;
#endif //_DEBUG
  sprintf(par->fnamePk,"default");
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

  par->do_lensing=0;
  par->oi_lens=NULL;
  
  par->do_gals=0;
  par->n_gals=-1;
  for(ii=0;ii<NPOP_MAX;ii++) {
    sprintf(par->fnameBzGals[ii],"default");
    sprintf(par->fnameNzGals[ii],"default");
    par->spline_gals_bz[ii]=NULL;
    par->spline_gals_nz[ii]=NULL;
    par->intacc_gals[ii]=NULL;
    par->shear_gals[ii]=0;
  }
  par->gals=NULL;
  par->nsources_this=NULL;

  return par;
}

static void conf_read_string(config_t *conf,char *secname,char *varname,char *out)
{
  int stat;
  char fullpath[256];
  const char *str;
  sprintf(fullpath,"%s.%s",secname,varname);
  stat=config_lookup_string(conf,fullpath,&str);
  if(stat==CONFIG_FALSE)
    report_error(1,"Couldn't read variable %s\n",fullpath);
  sprintf(out,"%s",str);
}

static void conf_read_double(config_t *conf,char *secname,char *varname,double *out)
{
  int stat;
  char fullpath[256];
  sprintf(fullpath,"%s.%s",secname,varname);
  stat=config_lookup_float(conf,fullpath,out);
  if(stat==CONFIG_FALSE)
    report_error(1,"Couldn't read variable %s\n",fullpath);
}

static void conf_read_int(config_t *conf,char *secname,char *varname,int *out)
{
  int stat;
  char fullpath[256];
  sprintf(fullpath,"%s.%s",secname,varname);
  stat=config_lookup_int(conf,fullpath,out);
  if(stat==CONFIG_FALSE)
    report_error(1,"Couldn't read variable %s\n",fullpath);
}

static void conf_read_bool(config_t *conf,char *secname,char *varname,int *out)
{
  int stat;
  char fullpath[256];
  sprintf(fullpath,"%s.%s",secname,varname);
  stat=config_lookup_bool(conf,fullpath,out);
  if(stat==CONFIG_FALSE)
    report_error(1,"Couldn't read variable %s\n",fullpath);
}

ParamCoLoRe *read_run_params(char *fname)
{
  int stat,ii,i_dum,found;
  char c_dum[256]="default";
  config_setting_t *cset;
  ParamCoLoRe *par=param_colore_new();
  config_t *conf=malloc(sizeof(config_t));
  config_init(conf);

  config_set_options(conf,CONFIG_OPTION_AUTOCONVERT);
  stat=config_read_file(conf,fname);
  if(stat==CONFIG_FALSE)
    error_open_file(fname);

  conf_read_string(conf,"global","prefix_out",par->prefixOut);
  conf_read_string(conf,"global","pk_filename",par->fnamePk);
  conf_read_double(conf,"global","r_smooth",&(par->r2_smooth));
  conf_read_double(conf,"global","z_min",&(par->z_min));
  conf_read_double(conf,"global","z_max",&(par->z_max));
  conf_read_int(conf,"global","n_grid",&(par->n_grid));
  conf_read_bool(conf,"global","output_density",&(par->output_density));
  conf_read_int(conf,"global","seed",&i_dum);
  par->seed_rng=i_dum;
  conf_read_string(conf,"global","output_format",c_dum);
  if(!strcmp(c_dum,"HDF5")) {
#ifdef _HAVE_HDF5
    par->output_format=2;
#else //_HAVE_HDF5
    report_error(1,"HDF5 format not supported\n");
#endif //_HAVE_HDF5
  }
  else if(!strcmp(c_dum,"FITS")) {
#ifdef _HAVE_FITS
    par->output_format=1;
#else //_HAVE_FITS
    report_error(1,"FITS format not supported\n");
#endif //_HAVE_FITS
  }
  else if(!strcmp(c_dum,"ASCII"))
    par->output_format=0;
  else
    report_error(1,"Unrecognized format %s\n",c_dum);

  conf_read_double(conf,"cosmo_par","omega_M",&(par->OmegaM));
  conf_read_double(conf,"cosmo_par","omega_L",&(par->OmegaL));
  conf_read_double(conf,"cosmo_par","omega_B",&(par->OmegaB));
  conf_read_double(conf,"cosmo_par","h",&(par->hhub));
  conf_read_double(conf,"cosmo_par","w",&(par->weos));
  conf_read_double(conf,"cosmo_par","ns",&(par->n_scal));
  conf_read_double(conf,"cosmo_par","sigma_8",&(par->sig8));

  //Get number of galaxy populations
  par->n_gals=0;
  found=1;
  while(found) {
    sprintf(c_dum,"gals%d",par->n_gals+1);
    cset=config_lookup(conf,c_dum);
    if(cset==NULL)
      found=0;
    else
      par->n_gals++;
  }
  if(par->n_gals>NPOP_MAX)
    report_error(1,"You're asking for too many populations %d! Enlarge NPOP_MAX\n",par->n_gals);
  for(ii=0;ii<par->n_gals;ii++) {
    sprintf(c_dum,"gals%d",ii+1);
    conf_read_string(conf,c_dum,"nz_filename",par->fnameNzGals[ii]);
    conf_read_string(conf,c_dum,"bias_filename",par->fnameBzGals[ii]);
    conf_read_bool(conf,c_dum,"include_shear",&(par->shear_gals[ii]));
    if(par->shear_gals[ii])
      par->do_lensing=1;
  }
  if(par->n_gals>0)
    par->do_gals=1;

  if(par->do_gals) {
    par->gals=my_malloc(par->n_gals*sizeof(Gal *));
    par->nsources_this=my_malloc(par->n_gals*sizeof(lint));
  }

#ifdef _DEBUG
  sprintf(c_dum,"%s_node%d.dbg",par->prefixOut,NodeThis);
  par->f_dbg=fopen(c_dum,"w");
  if(par->f_dbg==NULL) error_open_file(c_dum);
  if(NodeThis==0) {
    sprintf(c_dum,"%s_params.cfg",par->prefixOut);
    config_write_file(conf,c_dum);
  }
#endif //_DEBUG

  config_destroy(conf);

  if(par->r2_smooth>0) {
    par->r2_smooth=pow(par->r2_smooth,2);
    par->do_smoothing=1;
  }
  else
    par->do_smoothing=0;
  cosmo_set(par);
  init_fftw(par);
  if(par->do_lensing)
    par->oi_lens=alloc_onion_info(par,NSIDE_ONION_BASE,par->l_box/par->n_grid);

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
  if(par->do_gals)
    print_info("  %d galaxy populations\n",par->n_gals);
  if(par->do_lensing)
    print_info("  Will include lensing shear\n");
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
  int i_pop;
  char fname[256];
  
  if(par->do_gals) {
    if(NodeThis==0) timer(0);
    if(par->output_format==2) { //HDF5
#ifdef _HAVE_HDF5
      hid_t file_id,gal_types[5];
      hsize_t chunk_size=128;
      size_t dst_offset[4]={HOFFSET(Gal,ra),HOFFSET(Gal,dec),HOFFSET(Gal,z0),HOFFSET(Gal,dz_rsd)};
      const char *names[4]={"RA" ,"DEC","Z_COSMO","DZ_RSD"};
      char *tunit[4]=      {"DEG","DEG","NA"     ,"NA"    };
      gal_types[0]=H5T_NATIVE_FLOAT;
      gal_types[1]=H5T_NATIVE_FLOAT;
      gal_types[2]=H5T_NATIVE_FLOAT;
      gal_types[3]=H5T_NATIVE_FLOAT;
      
      print_info("*** Writing catalog (HDF5)\n");
      sprintf(fname,"%s_%d.h5",par->prefixOut,NodeThis);
      
      //Create file
      file_id=H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      //Write table for each galaxy type
      for(i_pop=0;i_pop<par->n_gals;i_pop++) {
	char table_title[256],table_name[256];
	sprintf(table_title,"sources%d_data",i_pop);
	sprintf(table_name,"/sources%d",i_pop);
	H5TBmake_table(table_title,file_id,table_name,4,par->nsources_this[i_pop],sizeof(Gal),
		       names,dst_offset,gal_types,chunk_size,NULL,0,par->gals[i_pop]);
	H5LTset_attribute_string(file_id,table_name,"FIELD_0_UNITS",tunit[0]);
	H5LTset_attribute_string(file_id,table_name,"FIELD_1_UNITS",tunit[1]);
	H5LTset_attribute_string(file_id,table_name,"FIELD_2_UNITS",tunit[2]);
	H5LTset_attribute_string(file_id,table_name,"FIELD_3_UNITS",tunit[3]);
      }
      
      //End file
      H5Fclose(file_id);
#else //_HAVE_HDF5
      printf("HDF5 not supported\n");
#endif //_HAVE_HDF5
    }
    else if(par->output_format==1) { //FITS
#ifdef _HAVE_FITS
      long ii,nrw=0;
      int status=0;
      fitsfile *fptr;
      int *type_arr;
      float *ra_arr,*dec_arr,*z0_arr,*rsd_arr;
      char *ttype[]={"RA" ,"DEC","Z_COSMO","DZ_RSD","TYPE"};
      char *tform[]={"1E" ,"1E" ,"1E"     ,"1E"    ,"1J"  };
      char *tunit[]={"DEG","DEG","NA"     ,"NA"    ,"NA"  };
      
      print_info("*** Writing catalog (FITS)\n");
      sprintf(fname,"!%s_%d.fits",par->prefixOut,NodeThis);
      
      fits_create_file(&fptr,fname,&status);
      fits_create_tbl(fptr,BINARY_TBL,0,5,ttype,tform,tunit,NULL,&status);
      
      fits_get_rowsize(fptr,&nrw,&status);
      ra_arr=my_malloc(nrw*sizeof(float));
      dec_arr=my_malloc(nrw*sizeof(float));
      z0_arr=my_malloc(nrw*sizeof(float));
      rsd_arr=my_malloc(nrw*sizeof(float));
      type_arr=my_malloc(nrw*sizeof(int));
      
      long offset=0;
      for(i_pop=0;i_pop<par->n_gals;i_pop++) {
	long row_here=0;
	while(row_here<par->nsources_this[i_pop]) {
	  long nrw_here;
	  if(row_here+nrw>par->nsources_this[i_pop])
	    nrw_here=par->nsources_this[i_pop]-row_here;
	  else
	    nrw_here=nrw;
	  
	  for(ii=0;ii<nrw_here;ii++) {
	    ra_arr[ii]=par->gals[i_pop][row_here+ii].ra;
	    dec_arr[ii]=par->gals[i_pop][row_here+ii].dec;
	    z0_arr[ii]=par->gals[i_pop][row_here+ii].z0;
	    rsd_arr[ii]=par->gals[i_pop][row_here+ii].dz_rsd;
	    type_arr[ii]=i_pop;
	  }
	  fits_write_col(fptr,TFLOAT,1,offset+row_here+1,1,nrw_here,ra_arr,&status);
	  fits_write_col(fptr,TFLOAT,2,offset+row_here+1,1,nrw_here,dec_arr,&status);
	  fits_write_col(fptr,TFLOAT,3,offset+row_here+1,1,nrw_here,z0_arr,&status);
	  fits_write_col(fptr,TFLOAT,4,offset+row_here+1,1,nrw_here,rsd_arr,&status);
	  fits_write_col(fptr,TINT  ,5,offset+row_here+1,1,nrw_here,type_arr,&status);
	  
	  row_here+=nrw_here;
	}
	offset+=par->nsources_this[i_pop];
      }
      fits_close_file(fptr,&status);
      free(ra_arr);
      free(dec_arr);
      free(z0_arr);
      free(rsd_arr);
      free(type_arr);
#else //_HAVE_FITS
      printf("FITS not supported\n");
#endif //_HAVE_FITS
    }
    else   {
      print_info("*** Writing catalog (ASCII)\n");
      sprintf(fname,"%s_%d.txt",par->prefixOut,NodeThis);
      
      lint jj;
      FILE *fil=fopen(fname,"w");
      if(fil==NULL) error_open_file(fname);
      fprintf(fil,"#[1]RA, [2]dec, [3]z0, [4]dz_RSD, [5]type\n");
      for(i_pop=0;i_pop<par->n_gals;i_pop++) {
	for(jj=0;jj<par->nsources_this[i_pop];jj++) {
	  fprintf(fil,"%E %E %E %E %d\n",
		  par->gals[i_pop][jj].ra,par->gals[i_pop][jj].dec,
		  par->gals[i_pop][jj].z0,par->gals[i_pop][jj].dz_rsd,
		  i_pop);
	}
      }
      fclose(fil);
    }
    if(NodeThis==0) timer(2);
    print_info("\n");
  }
}

void param_colore_free(ParamCoLoRe *par)
{
  int ii;

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
  if(par->do_lensing)
    free_onion_info(par->oi_lens);

  if(par->do_gals) {
    for(ii=0;ii<par->n_gals;ii++) {
      gsl_spline_free(par->spline_gals_bz[ii]);
      gsl_spline_free(par->spline_gals_nz[ii]);
      gsl_interp_accel_free(par->intacc_gals[ii]);
      if(par->gals[ii]!=NULL)
	free(par->gals[ii]);
    }
    if(par->gals!=NULL)
      free(par->gals);
    free(par->nsources_this);
  }

#ifdef _DEBUG
  fclose(par->f_dbg);
#endif //_DEBUG

  end_fftw();
}
