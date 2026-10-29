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
  ParamCoLoRe *par=my_malloc(sizeof(ParamCoLoRe));

#ifdef _DEBUG
  par->f_dbg=NULL;
#endif //_DEBUG

  //Cosmological parameters
  // Background
  par->OmegaM=0.3;
  par->OmegaL=0.7;
  par->OmegaB=0.05;
  par->hhub=0.7;
  par->weos=-1.;
  par->n_scal=0.96;
  par->sig8=0.83;
  par->prefac_lensing=-1;
  // Powerspectra
  sprintf(par->fnamePk,"default");
  par->numk=0;
  par->logkmax=1;
  par->logkmin=-3;
  par->idlogk=100;
  par->logkarr=NULL;
  par->pkarr=NULL;

  //Density parameters
  // Density methods
  par->output_density=0;
  par->r2_smooth=2.0;
  par->do_smoothing=1;
  par->smooth_potential=0;
  par->dens_type=DENS_TYPE_LGNR;
  par->lpt_interp_type=INTERP_CIC;
  par->lpt_buffer_fraction=0.2;
  par->output_lpt=0;
  par->seed_rng=1234;
  par->z0_norm=0;
  par->zf_norm=0;
  // Box parameters
  par->n_grid=512;
  par->l_box=-1;
  par->z_snap=-1;
  par->growth_d1=-1;
  par->growth_d2=-1;
  par->growth_dv=-1;
  par->growth_fz=-1;
  par->nz_here=512;
  par->iz0_here=0;
  par->nz_max=512;
  par->nz_all=NULL;
  par->iz0_all=NULL;
  // Density grids
  par->grid_dens_f=NULL;
  par->grid_dens=NULL;
  par->grid_npot_f=NULL;
  par->grid_npot=NULL;
  par->grid_eta_f=NULL;
  par->grid_eta=NULL;
  par->edge_planes=NULL;
  par->sigma2_gauss=0;

  //IO parameters
  sprintf(par->prefixOut,"default");
  par->output_format=0;

  //Tracers
  par->do_srcs=0;
  par->n_srcs=-1;
  par->bias=NULL;
  par->threshold=NULL;
  par->ndens=NULL;
  par->dens_norm=NULL;
  par->nsources_this=NULL;
  par->cats=NULL;

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

/*
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
*/

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

ParamCoLoRe *read_run_params(char *fname,int test_memory)
{
  int ii,stat,i_dum,found;
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
  conf_read_double(conf,"global","l_box", &(par->l_box));
  conf_read_double(conf,"global","z_snap", &(par->z_snap));
  conf_read_bool(conf,"global","output_density",&(par->output_density));
  conf_read_int(conf,"global","seed",&i_dum);
  conf_read_double(conf,"field_par","r_smooth",&(par->r2_smooth));
  conf_read_bool(conf,"field_par","smooth_potential",&(par->smooth_potential));
  conf_read_int(conf,"field_par","n_grid",&(par->n_grid));
  conf_read_int(conf,"field_par","dens_type",&(par->dens_type));
  conf_read_double(conf,"field_par","lpt_buffer_fraction",&(par->lpt_buffer_fraction));
  conf_read_int(conf,"field_par","lpt_interp_type",&(par->lpt_interp_type));
  conf_read_int(conf,"field_par","output_lpt",&(par->output_lpt));

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
  if(par->n_srcs>0) {
    par->do_srcs=1;
    par->bias=my_malloc(par->n_srcs*sizeof(double));
    par->threshold=my_malloc(par->n_srcs*sizeof(double));
    par->ndens=my_malloc(par->n_srcs*sizeof(double));
    par->dens_norm=my_malloc(par->n_srcs*sizeof(double));
    par->nsources_this=my_calloc(par->n_srcs,sizeof(long));
    par->cats=my_malloc(par->n_srcs*sizeof(CatalogCartesian *));
    for(ii=0;ii<par->n_srcs;ii++) {
      sprintf(c_dum,"srcs%d",ii+1);
      conf_read_double(conf,c_dum,"bias",&(par->bias[ii]));
      conf_read_double(conf,c_dum,"threshold",&(par->threshold[ii]));
      conf_read_double(conf,c_dum,"ndens",&(par->ndens[ii]));
    }
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

  init_fftw(par);

  cosmo_set(par);

  get_max_memory(par,test_memory);

  print_info("\n");

  double dk=2*M_PI/par->l_box;
  print_info("Run parameters: \n");
  print_info("  L_box = %.3lf Mpc/h, N_grid = %d \n",(double)(par->l_box),par->n_grid);
  print_info("  Scales resolved: %.3lE < k < %.3lE h/Mpc\n",dk,0.5*(par->n_grid-1)*dk);
  print_info("  Fourier-space resolution: dk = %.3lE h/Mpc\n",dk);
  print_info("  Real-space resolution: dx = %.3lE Mpc/h\n",par->l_box/par->n_grid);
  if(par->do_smoothing)
    print_info("  Density field pre-smoothed on scales: x_s = %.3lE Mpc/h\n",sqrt(par->r2_smooth));
  else
    print_info("  No extra smoothing\n");
  if(par->do_srcs)
    print_info("  %d galaxy populations\n",par->n_srcs);
  print_info("\n");

  if(test_memory) {
#ifdef _DEBUG
    fclose(par->f_dbg);
#endif //_DEBUG
    free(par);
    return NULL;
  }

  allocate_fftw(par);

  free(conf);

  return par;
}

void write_density_grid(ParamCoLoRe *par,char *prefix_dens)
{
  FILE *fo;
  char fname[256];
  int iz;
  int ngx=2*(par->n_grid/2+1);
  int size_flouble=sizeof(flouble);
  double lb=par->l_box;

  if(NodeThis==0) timer(0);
  print_info("*** Writing density field (native format)\n");
  sprintf(fname,"%s_dens_%s_%d.dat",par->prefixOut,prefix_dens,NodeThis);
  fo=fopen(fname,"wb");
  if(fo==NULL) error_open_file(fname);
  my_fwrite(&NNodes,sizeof(int),1,fo);
  my_fwrite(&size_flouble,sizeof(int),1,fo);
  my_fwrite(&lb,sizeof(double),1,fo);
  my_fwrite(&(par->n_grid),sizeof(int),1,fo);
  my_fwrite(&(par->nz_here),sizeof(int),1,fo);
  my_fwrite(&(par->iz0_here),sizeof(int),1,fo);
  for(iz=0;iz<par->nz_here;iz++) {
    int iy;
    for(iy=0;iy<par->n_grid;iy++) {
      long index0=ngx*((long)(iy+iz*par->n_grid));
      my_fwrite(&(par->grid_dens[index0]),sizeof(flouble),par->n_grid,fo);
    }
  }
  fclose(fo);
  if(NodeThis==0) timer(2);
  print_info("\n");
}

void write_eta_grid(ParamCoLoRe *par)
{
  FILE *fo;
  char fname[256];
  int iz;
  int ngx=2*(par->n_grid/2+1);
  int size_flouble=sizeof(flouble);
  double lb=par->l_box;

  if(NodeThis==0) timer(0);
  print_info("*** Writing density field (native format)\n");
  sprintf(fname,"%s_eta_%d.dat",par->prefixOut,NodeThis);
  fo=fopen(fname,"wb");
  if(fo==NULL) error_open_file(fname);
  my_fwrite(&NNodes,sizeof(int),1,fo);
  my_fwrite(&size_flouble,sizeof(int),1,fo);
  my_fwrite(&lb,sizeof(double),1,fo);
  my_fwrite(&(par->n_grid),sizeof(int),1,fo);
  my_fwrite(&(par->nz_here),sizeof(int),1,fo);
  my_fwrite(&(par->iz0_here),sizeof(int),1,fo);
  for(iz=0;iz<par->nz_here;iz++) {
    int iy;
    for(iy=0;iy<par->n_grid;iy++) {
      long index0=ngx*((long)(iy+iz*par->n_grid));
      my_fwrite(&(par->grid_eta[index0]),sizeof(flouble),par->n_grid,fo);
    }
  }
  fclose(fo);
  if(NodeThis==0) timer(2);
  print_info("\n");
}

typedef struct {
  int    np[6];
  double mass[6];
  double time;
  double redshift;
  int    flag_sfr;
  int    flag_feedback;
  unsigned int np_total[6];
  int    flag_cooling;
  int    num_files;
  double boxsize;
  double omega0;
  double omega_lambda;
  double hubble_param;
  int flag_stellarage;
  int flag_metals;
  unsigned int np_total_highword[6];
  int  flag_entropy_instead_u;
  int flag_gadgetformat;
  char fill[56];
} GadgetHeader;

void write_lpt(ParamCoLoRe *par,unsigned long long npart,flouble *x,flouble *y,flouble *z)
{
  GadgetHeader header;
  FILE *fo;
  char fname[256];
  unsigned long long ipart,np_total;
  unsigned long long np_total_expected=par->n_grid*((long)(par->n_grid*par->n_grid));

  sprintf(fname,"%s_lpt_out.%d",par->prefixOut,NodeThis);
  fo=fopen(fname,"w");
  if(fo==NULL) error_open_file(fname);

#ifdef _HAVE_MPI
  unsigned long long np_send=npart;
  MPI_Reduce(&np_send,&np_total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Bcast(&np_total,1,MPI_UNSIGNED_LONG_LONG,0,MPI_COMM_WORLD);
#else //_HAVE_MPI
  np_total=npart;
#endif //_HAVE_MPI

  if(np_total!=np_total_expected)
    report_error(1,"Only %llu particles found, but there should be %ull\n",np_total,np_total_expected);

  double m=27.7455*par->OmegaM*pow(par->l_box,3.)/np_total;
  memset(&header,0,sizeof(GadgetHeader));

  header.np[1]=npart;
  header.mass[1]=m;
  header.time=1.;
  header.redshift=0.;
  header.np_total[1]=(unsigned int)np_total;
  header.np_total_highword[1]=(unsigned int)(np_total >> 32);
  header.num_files=NNodes;
  header.boxsize=par->l_box;
  header.omega0=par->OmegaM;
  header.omega_lambda=par->OmegaL;
  header.hubble_param=par->hhub;
  header.flag_gadgetformat=1;

  int blklen=sizeof(GadgetHeader);
  my_fwrite(&blklen,sizeof(blklen),1,fo);
  my_fwrite(&header,sizeof(GadgetHeader),1,fo);
  my_fwrite(&blklen,sizeof(blklen),1,fo);

  float x0[3];
  // position
  blklen=npart*sizeof(float)*3;
  my_fwrite(&blklen,sizeof(blklen),1,fo);
  for(ipart=0;ipart<npart;ipart++) {
    x0[0]=x[ipart];
    x0[1]=y[ipart];
    x0[2]=z[ipart];
    my_fwrite(x0,sizeof(float),3,fo);
  }
  my_fwrite(&blklen,sizeof(blklen),1,fo);

  // velocity
  x0[0]=0; x0[1]=0; x0[2]=0;
  blklen=npart*sizeof(float)*3;
  my_fwrite(&blklen,sizeof(blklen),1,fo);
  for(ipart=0;ipart<npart;ipart++) {
    my_fwrite(x0,sizeof(float),3,fo);
  }
  my_fwrite(&blklen,sizeof(blklen),1,fo);

  // id
  blklen=npart*sizeof(unsigned long long);
  my_fwrite(&blklen,sizeof(blklen),1,fo);
  long long id0=(long long)(par->iz0_here*((long)(par->n_grid*par->n_grid)));
  for(ipart=0;ipart<npart;ipart++) {
    unsigned long long id_out=id0+ipart;
    my_fwrite(&id_out,sizeof(unsigned long long),1,fo);
  }
  my_fwrite(&blklen,sizeof(blklen),1,fo);

  fclose(fo);
}

static void write_catalog(ParamCoLoRe *par,int ipop)
{
  char fname[256];

  if(NodeThis==0) timer(0);
#ifdef _HAVE_FITS
  long ii,nrw=0;
  int status=0;
  fitsfile *fptr;
  int *type_arr;
  float *x_arr,*y_arr,*z_arr,*rsd_arr;
  int nfields=5;
  char *ttype[]={"TYPE","X"    ,"Y"    ,"Z"    ,"DX_RSD"};
  char *tform[]={"1J"  ,"1E"   ,"1E"   ,"1E"   ,"1E"};
  char *tunit[]={"NA"  ,"Mpc/h","Mpc/h","Mpc/h","Mpc/h"};

  print_info(" %d-th population (FITS)\n",ipop);
  sprintf(fname,"!%s_srcs_s%d_%d.fits",par->prefixOut,ipop+1,NodeThis);

  fits_create_file(&fptr,fname,&status);
  fits_create_tbl(fptr,BINARY_TBL,0,nfields,ttype,tform,tunit,NULL,&status);
  fits_update_key(fptr,TSTRING,"CONTENTS","Source catalog",NULL,&status);

  fits_get_rowsize(fptr,&nrw,&status);
  type_arr=my_malloc(nrw*sizeof(int));
  x_arr=my_malloc(nrw*sizeof(float));
  y_arr=my_malloc(nrw*sizeof(float));
  z_arr=my_malloc(nrw*sizeof(float));
  rsd_arr=my_malloc(nrw*sizeof(float));

  long row_here=0;
  while(row_here<par->nsources_this[ipop]) {
    long nrw_here;
    if(row_here+nrw>par->nsources_this[ipop])
      nrw_here=par->nsources_this[ipop]-row_here;
    else
      nrw_here=nrw;

    for(ii=0;ii<nrw_here;ii++) {
      type_arr[ii]=ipop;
      x_arr[ii]=par->cats[ipop]->pos[(row_here+ii)*NPOS_CC+0];
      y_arr[ii]=par->cats[ipop]->pos[(row_here+ii)*NPOS_CC+1];
      z_arr[ii]=par->cats[ipop]->pos[(row_here+ii)*NPOS_CC+2];
      rsd_arr[ii]=par->cats[ipop]->pos[(row_here+ii)*NPOS_CC+3];
    }
    fits_write_col(fptr,TINT  ,1,row_here+1,1,nrw_here,type_arr,&status);
    fits_write_col(fptr,TFLOAT,2,row_here+1,1,nrw_here,x_arr,&status);
    fits_write_col(fptr,TFLOAT,3,row_here+1,1,nrw_here,y_arr,&status);
    fits_write_col(fptr,TFLOAT,4,row_here+1,1,nrw_here,z_arr,&status);
    fits_write_col(fptr,TFLOAT,5,row_here+1,1,nrw_here,rsd_arr,&status);
    row_here+=nrw_here;
  }
  fits_close_file(fptr,&status);
  free(x_arr);
  free(y_arr);
  free(z_arr);
  free(rsd_arr);
#else //_HAVE_FITS
  printf("FITS not supported\n");
#endif //_HAVE_FITS
  if(NodeThis==0) timer(2);
}

void write_srcs(ParamCoLoRe *par)
{
  print_info("*** Writing source catalogs\n");
  int ipop;
  for(ipop=0;ipop<par->n_srcs;ipop++)
    write_catalog(par,ipop);
  print_info("\n");
}

void param_colore_free(ParamCoLoRe *par)
{
    int ii;
  free(par->nz_all);
  free(par->iz0_all);
  free(par->logkarr);
  free(par->pkarr);
  end_fftw(par);

  if(par->do_srcs) {
    free(par->bias);
    free(par->threshold);
    free(par->ndens);
    free(par->dens_norm);
    if(par->cats!=NULL) {
      for(ii=0;ii<par->n_srcs;ii++)
	catalog_cartesian_free(par->cats[ii]);
      free(par->cats);
    }
    free(par->nsources_this);
  }

#ifdef _DEBUG
  fclose(par->f_dbg);
#endif //_DEBUG
  free(par);
}
