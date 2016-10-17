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
  par->smooth_potential=0;
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
  par->grid_npot_f=NULL;
  par->grid_npot=NULL;
  par->grid_dumm=NULL;
  par->sigma2_gauss=0;

  par->do_lensing=0;
  par->nside_base=-1;
  par->n_beams_here=-1;
  par->oi_slices=NULL;
  par->oi_beams=NULL;
  par->oi_sl_dum=NULL;
  par->dens_slices=NULL;
  par->vrad_slices=NULL;
  par->p_xx_slices=NULL;
  par->p_xy_slices=NULL;
  par->p_yy_slices=NULL;
  par->dens_beams=NULL;
  par->vrad_beams=NULL;
  par->p_xx_beams=NULL;
  par->p_xy_beams=NULL;
  par->p_yy_beams=NULL;
  par->nsrc_beams=NULL;

  par->do_sources=0;
  par->do_imap=0;
  par->do_kappa=0;
  par->do_pred=0;
  par->n_srcs=-1;
  par->n_imap=-1;
  par->n_kappa=-1;
  par->nside_kappa=-1;
  for(ii=0;ii<NPOP_MAX;ii++) {
    sprintf(par->fnameBzSrcs[ii],"default");
    sprintf(par->fnameNzSrcs[ii],"default");
    par->spline_srcs_bz[ii]=NULL;
    par->spline_srcs_nz[ii]=NULL;
    par->intacc_srcs[ii]=NULL;
    par->shear_srcs[ii]=0;
    sprintf(par->fnameBzImap[ii],"default");
    sprintf(par->fnameTzImap[ii],"default");
    sprintf(par->fnameNuImap[ii],"default");
    par->spline_imap_bz[ii]=NULL;
    par->spline_imap_tz[ii]=NULL;
    par->intacc_imap[ii]=NULL;
    par->nside_imap[ii]=-1;
    par->z_kappa_out[ii]=-1;
    par->nu0_imap[ii]=-1;
  }
  par->srcs=NULL;
  par->imap=NULL;
  par->kmap=NULL;
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

static void conf_read_double_array(config_t *conf,char *secname,char *varname,
				   double *out,int *nel,int nmax)
{
  int n_elem,index;
  config_setting_t *arr;
  char fullpath[256];
  sprintf(fullpath,"%s.%s",secname,varname);
  arr=config_lookup(conf,fullpath);
  if(arr==NULL)
    report_error(1,"Couldn't read variable %s\n",fullpath);
  n_elem=config_setting_length(arr);
  if(n_elem==0)
    report_error(1,"Couldn't read variable %s\n",fullpath);
  if(n_elem>nmax)
    report_error(1,"Too many elements in %s (%d > %d)\n",fullpath,n_elem,nmax);

  *nel=n_elem;
  for(index=0;index<n_elem;index++) {
    out[index]=config_setting_get_float_elem(arr,index);
  }
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
  conf_read_bool(conf,"global","smooth_potential",&(par->smooth_potential));
  conf_read_double(conf,"global","z_min",&(par->z_min));
  conf_read_double(conf,"global","z_max",&(par->z_max));
  conf_read_int(conf,"global","n_grid",&(par->n_grid));
  conf_read_bool(conf,"global","output_density",&(par->output_density));
  conf_read_int(conf,"global","seed",&i_dum);
  conf_read_bool(conf,"global","write_pred",&(par->do_pred));
  if (par->do_pred) 
    conf_read_double(conf,"global","pred_dz",&(par->pred_dz));

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
  par->n_srcs=0;
  found=1;
  while(found) {
    sprintf(c_dum,"srcs%d",par->n_srcs+1);
    cset=config_lookup(conf,c_dum);
    if(cset==NULL)
      found=0;
    else
      par->n_srcs++;
  }
  if(par->n_srcs>NPOP_MAX)
    report_error(1,"You're asking for too many populations %d! Enlarge NPOP_MAX\n",par->n_srcs);
  for(ii=0;ii<par->n_srcs;ii++) {
    sprintf(c_dum,"srcs%d",ii+1);
    conf_read_string(conf,c_dum,"nz_filename",par->fnameNzSrcs[ii]);
    conf_read_string(conf,c_dum,"bias_filename",par->fnameBzSrcs[ii]);
    conf_read_bool(conf,c_dum,"include_shear",&(par->shear_srcs[ii]));
    if(par->shear_srcs[ii])
      par->do_lensing=1;
  }
  if(par->n_srcs>0)
    par->do_sources=1;

  //Get number of intensity mapping species
  par->n_imap=0;
  found=1;
  while(found) {
    sprintf(c_dum,"imap%d",par->n_imap+1);
    cset=config_lookup(conf,c_dum);
    if(cset==NULL)
      found=0;
    else
      par->n_imap++;
  }
  if(par->n_imap>NPOP_MAX)
    report_error(1,"You're asking for too many populations %d! Enlarge NPOP_MAX\n",par->n_imap);
  for(ii=0;ii<par->n_imap;ii++) {
    sprintf(c_dum,"imap%d",ii+1);
    conf_read_string(conf,c_dum,"tbak_filename",par->fnameTzImap[ii]);
    conf_read_string(conf,c_dum,"bias_filename",par->fnameBzImap[ii]);
    conf_read_string(conf,c_dum,"freq_list",par->fnameNuImap[ii]);
    conf_read_double(conf,c_dum,"freq_rest",&(par->nu0_imap[ii]));
    conf_read_int(conf,c_dum,"nside",&(par->nside_imap[ii]));
  }
  if(par->n_imap>0)
    par->do_imap=1;

  //Kappa maps
  cset=config_lookup(conf,"kappa");
  if(cset!=NULL) {
    par->do_kappa=1;
    par->do_lensing=1;
    conf_read_double_array(conf,"kappa","z_out",par->z_kappa_out,&(par->n_kappa),NPOP_MAX);
    conf_read_int(conf,"kappa","nside",&(par->nside_kappa));
  }

  if(par->do_sources) {
    par->srcs=my_malloc(par->n_srcs*sizeof(Src *));
    par->nsources_this=my_malloc(par->n_srcs*sizeof(lint));
  }

  if(par->do_imap) {
    par->imap=my_malloc(par->n_imap*sizeof(HealpixShells *));
    for(ii=0;ii<par->n_imap;ii++) {
      FILE *fnu=fopen(par->fnameNuImap[ii],"r");
      if(fnu==NULL) error_open_file(par->fnameNuImap[ii]);
      par->imap[ii]=new_hp_shell(par->nside_imap[ii],linecount(fnu));
      fclose(fnu);
    }
  }

  if(par->do_kappa)
    par->kmap=new_hp_shell(par->nside_kappa,par->n_kappa);

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

  par->nside_base=1;
  while(2*par->nside_base*par->nside_base<NNodes)
    par->nside_base*=2;

  par->oi_slices=alloc_onion_info_slices(par);
  par->oi_sl_dum=alloc_onion_info_slices(par);
  par->oi_beams=alloc_onion_info_beams(par);
  par->nside_base=par->oi_slices->nside_arr[0];
  alloc_slices(par);

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
  if(par->do_sources)
    print_info("  %d galaxy populations\n",par->n_srcs);
  if(par->do_imap)
    print_info("  %d intensity mapping species\n",par->n_imap);
  if(par->do_kappa)
    print_info("  %d lensing source planes\n",par->n_kappa);
  if(par->do_lensing)
    print_info("  Will include lensing shear\n");
  print_info("\n");
  
  return par;
}

void write_grids(ParamCoLoRe *par)
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

void write_imap(ParamCoLoRe *par)
{
  int i_pop;
  char fname[256];

  if(par->do_imap) {
    if(NodeThis==0) timer(0);
    print_info("*** Writing intensity maps\n");
    for(i_pop=0;i_pop<par->n_imap;i_pop++) {
      int ir;
      long npx=he_nside2npix(par->imap[i_pop]->nside);
      flouble *map_write=my_malloc(npx*sizeof(flouble));
      for(ir=0;ir<par->imap[i_pop]->nr;ir++) {
	long ip;
	long ir_t=ir*par->imap[i_pop]->num_pix;

	//Write local pixels to dummy map
	memset(map_write,0,npx*sizeof(flouble));
	sprintf(fname,"!%s_imap_s%d_nu%03d.fits",par->prefixOut,i_pop,ir);
	for(ip=0;ip<npx;ip++) {
	  int id_pix=par->imap[i_pop]->listpix[ip];
	  if(id_pix>0)
	    map_write[ip]+=par->imap[i_pop]->data[ir_t+id_pix];
	}

	//Collect all dummy maps
#ifdef _HAVE_MPI
	if(NodeThis==0)
	  MPI_Reduce(MPI_IN_PLACE,map_write,npx,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
	else
	  MPI_Reduce(map_write   ,NULL     ,npx,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
#endif //_HAVE_MPI

	//Write dummy map
	if(NodeThis==0)
	  he_write_healpix_map(&map_write,1,par->imap[i_pop]->nside,fname);
      }
      free(map_write);
    }
    if(NodeThis==0) timer(2);
    print_info("\n");
  }
}

void write_kappa(ParamCoLoRe *par)
{
  if(par->do_kappa) {
    int ir;
    char fname[256];
    long npx=he_nside2npix(par->nside_kappa);
    flouble *map_write=my_malloc(npx*sizeof(flouble));
    if(NodeThis==0) timer(0);
    print_info("*** Writing kappa source maps\n");
    for(ir=0;ir<par->kmap->nr;ir++) {
      long ip;
      long ir_t=ir*par->kmap->num_pix;
      
      //Write local pixels to dummy map
      memset(map_write,0,npx*sizeof(flouble));
      sprintf(fname,"!%s_kappa_z%03d.fits",par->prefixOut,ir);
      for(ip=0;ip<npx;ip++) {
	int id_pix=par->kmap->listpix[ip];
	if(id_pix>0)
	  map_write[ip]+=par->kmap->data[ir_t+id_pix];
      }

      //Collect all dummy maps
#ifdef _HAVE_MPI
      if(NodeThis==0)
	MPI_Reduce(MPI_IN_PLACE,map_write,npx,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
      else
	MPI_Reduce(map_write   ,NULL     ,npx,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
#endif //_HAVE_MPI
      
      //Write dummy map
      if(NodeThis==0)
	he_write_healpix_map(&map_write,1,par->nside_kappa,fname);
    }
    free(map_write);
    if(NodeThis==0) timer(2);
    print_info("\n");
  }
}

void write_catalog(ParamCoLoRe *par)
{
  int i_pop;
  char fname[256];

  if(par->do_sources) {
    if(NodeThis==0) timer(0);
    if(par->output_format==2) { //HDF5
#ifdef _HAVE_HDF5
      hid_t file_id,gal_types[6];
      hsize_t chunk_size=128;
      size_t dst_offset[6]={HOFFSET(Src,ra),HOFFSET(Src,dec),HOFFSET(Src,z0),HOFFSET(Src,dz_rsd),HOFFSET(Src,e1),HOFFSET(Src,e2)};
      const char *names[6]={"RA" ,"DEC","Z_COSMO","DZ_RSD","E1","E2"};
      char *tunit[6]=      {"DEG","DEG","NA"     ,"NA"    ,"NA","NA"};
      gal_types[0]=H5T_NATIVE_FLOAT;
      gal_types[1]=H5T_NATIVE_FLOAT;
      gal_types[2]=H5T_NATIVE_FLOAT;
      gal_types[3]=H5T_NATIVE_FLOAT;
      gal_types[4]=H5T_NATIVE_FLOAT;
      gal_types[5]=H5T_NATIVE_FLOAT;

      print_info("*** Writing catalog (HDF5)\n");
      sprintf(fname,"%s_srcs_%d.h5",par->prefixOut,NodeThis);

      //Create file
      file_id=H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      //Write table for each galaxy type
      for(i_pop=0;i_pop<par->n_srcs;i_pop++) {
	char table_title[256],table_name[256];
	sprintf(table_title,"sources%d_data",i_pop);
	sprintf(table_name,"/sources%d",i_pop);
	H5TBmake_table(table_title,file_id,table_name,6,par->nsources_this[i_pop],sizeof(Src),
		       names,dst_offset,gal_types,chunk_size,NULL,0,par->srcs[i_pop]);
	H5LTset_attribute_string(file_id,table_name,"FIELD_0_UNITS",tunit[0]);
	H5LTset_attribute_string(file_id,table_name,"FIELD_1_UNITS",tunit[1]);
	H5LTset_attribute_string(file_id,table_name,"FIELD_2_UNITS",tunit[2]);
	H5LTset_attribute_string(file_id,table_name,"FIELD_3_UNITS",tunit[3]);
	H5LTset_attribute_string(file_id,table_name,"FIELD_4_UNITS",tunit[4]);
	H5LTset_attribute_string(file_id,table_name,"FIELD_5_UNITS",tunit[5]);
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
      float *ra_arr,*dec_arr,*z0_arr,*rsd_arr,*e1_arr,*e2_arr;
      int nfields=5;
      char *ttype[]={"TYPE","RA" ,"DEC","Z_COSMO","DZ_RSD","E1","E2"};
      char *tform[]={"1J"  ,"1E" ,"1E" ,"1E"     ,"1E"    ,"1E","1E"};
      char *tunit[]={"NA"  ,"DEG","DEG","NA"     ,"NA"    ,"NA","NA"};
      if(par->do_lensing)
	nfields=7;

      print_info("*** Writing catalog (FITS)\n");
      sprintf(fname,"!%s_srcs_%d.fits",par->prefixOut,NodeThis);

      fits_create_file(&fptr,fname,&status);
      fits_create_tbl(fptr,BINARY_TBL,0,nfields,ttype,tform,tunit,NULL,&status);

      fits_get_rowsize(fptr,&nrw,&status);
      type_arr=my_malloc(nrw*sizeof(int));
      ra_arr=my_malloc(nrw*sizeof(float));
      dec_arr=my_malloc(nrw*sizeof(float));
      z0_arr=my_malloc(nrw*sizeof(float));
      rsd_arr=my_malloc(nrw*sizeof(float));
      e1_arr=my_malloc(nrw*sizeof(float));
      e2_arr=my_malloc(nrw*sizeof(float));

      long offset=0;
      for(i_pop=0;i_pop<par->n_srcs;i_pop++) {
	long row_here=0;
	while(row_here<par->nsources_this[i_pop]) {
	  long nrw_here;
	  if(row_here+nrw>par->nsources_this[i_pop])
	    nrw_here=par->nsources_this[i_pop]-row_here;
	  else
	    nrw_here=nrw;

	  for(ii=0;ii<nrw_here;ii++) {
	    type_arr[ii]=i_pop;
	    ra_arr[ii]=par->srcs[i_pop][row_here+ii].ra;
	    dec_arr[ii]=par->srcs[i_pop][row_here+ii].dec;
	    z0_arr[ii]=par->srcs[i_pop][row_here+ii].z0;
	    rsd_arr[ii]=par->srcs[i_pop][row_here+ii].dz_rsd;
	    if(par->do_lensing) {
	      e1_arr[ii]=par->srcs[i_pop][row_here+ii].e1;
	      e2_arr[ii]=par->srcs[i_pop][row_here+ii].e2;
	    }
	  }
	  fits_write_col(fptr,TINT  ,1,offset+row_here+1,1,nrw_here,type_arr,&status);
	  fits_write_col(fptr,TFLOAT,2,offset+row_here+1,1,nrw_here,ra_arr,&status);
	  fits_write_col(fptr,TFLOAT,3,offset+row_here+1,1,nrw_here,dec_arr,&status);
	  fits_write_col(fptr,TFLOAT,4,offset+row_here+1,1,nrw_here,z0_arr,&status);
	  fits_write_col(fptr,TFLOAT,5,offset+row_here+1,1,nrw_here,rsd_arr,&status);
	  if(par->do_lensing) {
	    fits_write_col(fptr,TFLOAT,6,offset+row_here+1,1,nrw_here,e1_arr,&status);
	    fits_write_col(fptr,TFLOAT,7,offset+row_here+1,1,nrw_here,e2_arr,&status);
	  }

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
      free(e1_arr);
      free(e2_arr);
#else //_HAVE_FITS
      printf("FITS not supported\n");
#endif //_HAVE_FITS
    }
    else   {
      print_info("*** Writing catalog (ASCII)\n");
      sprintf(fname,"%s_srcs_%d.txt",par->prefixOut,NodeThis);

      lint jj;
      FILE *fil=fopen(fname,"w");
      if(fil==NULL) error_open_file(fname);
      fprintf(fil,"#[1]type [2]RA, [3]dec, [4]z0, [5]dz_RSD ");
      if(par->do_lensing)
	fprintf(fil,"#[6]e1, [7]e2\n");
      else
	fprintf(fil,"\n");
      for(i_pop=0;i_pop<par->n_srcs;i_pop++) {
	for(jj=0;jj<par->nsources_this[i_pop];jj++) {
	  fprintf(fil,"%d %E %E %E %E ",
		  i_pop,par->srcs[i_pop][jj].ra,par->srcs[i_pop][jj].dec,
		  par->srcs[i_pop][jj].z0,par->srcs[i_pop][jj].dz_rsd);
	  if(par->do_lensing)
	    fprintf(fil,"%E %E \n",par->srcs[i_pop][jj].e1,par->srcs[i_pop][jj].e2);
	  else
	    fprintf(fil,"\n");
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

  free_beams(par);
  free_onion_info(par->oi_slices);
  free_onion_info(par->oi_sl_dum);
  for(ii=0;ii<par->n_beams_here;ii++)
    free_onion_info(par->oi_beams[ii]);
  free(par->oi_beams);

  if(par->do_sources) {
    for(ii=0;ii<par->n_srcs;ii++) {
      gsl_spline_free(par->spline_srcs_bz[ii]);
      gsl_spline_free(par->spline_srcs_nz[ii]);
      gsl_interp_accel_free(par->intacc_srcs[ii]);
      if(par->srcs[ii]!=NULL)
      	free(par->srcs[ii]);
    }
    if(par->srcs!=NULL)
      free(par->srcs);
    free(par->nsources_this);
  }

  if(par->do_imap) {
    for(ii=0;ii<par->n_imap;ii++) {
      gsl_spline_free(par->spline_imap_bz[ii]);
      gsl_spline_free(par->spline_imap_tz[ii]);
      gsl_interp_accel_free(par->intacc_imap[ii]);
      free_hp_shell(par->imap[ii]);
    }
    if(par->imap!=NULL)
      free(par->imap);
  }

  if(par->do_kappa) {
    free_hp_shell(par->kmap);
  }

#ifdef _DEBUG
  fclose(par->f_dbg);
#endif //_DEBUG
}
