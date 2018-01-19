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

  //Cosmological parameters
  // Background
  par->OmegaM=0.3;
  par->OmegaL=0.7;
  par->OmegaB=0.05;
  par->hhub=0.7;
  par->weos=-1.;
  par->n_scal=0.96;
  par->sig8=0.83;
  par->fgrowth_0=-1;
  par->hubble_0=-1;
  par->prefac_lensing=-1;
  par->z_max=1.5;
  par->z_min=0.5;
  par->r_max=-1;
  par->r_min=-1;
  par->glob_idr=-1;
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
  par->sigma2_gauss=0;

  //IO parameters
  sprintf(par->prefixOut,"default");
  par->output_format=0;
  par->do_pred=0;
  par->do_pred=1;

  //Tracers
  par->do_srcs=0;
  par->do_skewers=0;
  par->do_lensing=0;
  par->do_isw=0;
  par->do_imap=0;
  par->do_kappa=0;
  par->do_isw=0;
  par->n_srcs=-1;
  par->n_imap=-1;
  par->n_kappa=-1;
  par->n_isw=-1;
  par->nside_kappa=-1;
  par->nside_isw=-1;
  for(ii=0;ii<NPOP_MAX;ii++) {
    sprintf(par->fnameBzSrcs[ii],"default");
    sprintf(par->fnameNzSrcs[ii],"default");
    par->srcs_bz_arr[ii]=NULL;
    par->srcs_nz_arr[ii]=NULL;
    par->srcs_norm_arr[ii]=NULL;
    par->shear_srcs[ii]=0;
    par->skw_srcs[ii]=0;
    sprintf(par->fnameBzImap[ii],"default");
    sprintf(par->fnameTzImap[ii],"default");
    sprintf(par->fnameNuImap[ii],"default");
    par->imap_bz_arr[ii]=NULL;
    par->imap_tz_arr[ii]=NULL;
    par->imap_norm_arr[ii]=NULL;
    par->nside_imap[ii]=-1;
    par->nu0_imap[ii]=-1;
  }
  for(ii=0;ii<NPLANES_MAX;ii++) {
    par->z_kappa_out[ii]=-1;
    par->z_isw_out[ii]=-1;
  }
  par->cats_c=NULL;
  par->nsources_c_this=NULL;
  par->cats=NULL;
  par->nsources_this=NULL;
  par->imap=NULL;
  par->kmap=NULL;
  par->pd_map=NULL;

  //Beam distribution
  par->nside_base=-1;
  par->npix_base=-1;
  par->need_beaming=0;

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

static int choose_nside_base(void)
{
  int nside_base=2;
  int enough=0;

  while(!enough) {
    long npix=he_nside2npix(nside_base);
    if(npix%NNodes==0)
      enough=1;
    else {
      int pernode=npix/NNodes;
      double load_diff=(pernode+1.)/pernode;
      if(load_diff<1.2)
	enough=1;
      else
	nside_base*=2;
    }
  }

  return nside_base;
}

ParamCoLoRe *read_run_params(char *fname,int test_memory)
{
  int stat,ii,i_dum,found;
  char c_dum[256]="default";
  config_setting_t *cset;
  ParamCoLoRe *par=param_colore_new();
  config_t *conf=malloc(sizeof(config_t));

  par->nside_base=choose_nside_base();
  if(par->nside_base>NSIDE_MAX_HPX)
    report_error(1,"Can't go beyond nside=%d\n",NSIDE_MAX_HPX);
  par->npix_base=he_nside2npix(par->nside_base);

  config_init(conf);

  config_set_options(conf,CONFIG_OPTION_AUTOCONVERT);
  stat=config_read_file(conf,fname);
  if(stat==CONFIG_FALSE)
    error_open_file(fname);

  conf_read_string(conf,"global","prefix_out",par->prefixOut);
  conf_read_string(conf,"global","pk_filename",par->fnamePk);
  conf_read_double(conf,"global","z_min",&(par->z_min));
  conf_read_double(conf,"global","z_max",&(par->z_max));
  conf_read_bool(conf,"global","output_density",&(par->output_density));
  conf_read_int(conf,"global","seed",&i_dum);
  conf_read_bool(conf,"global","write_pred",&(par->do_pred));
  if (par->do_pred) {
    conf_read_bool(conf,"global","just_write_pred",&(par->just_do_pred));
    conf_read_double(conf,"global","pred_dz",&(par->pred_dz));
  }
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
  for(ii=0;ii<par->n_srcs;ii++) {
    sprintf(c_dum,"srcs%d",ii+1);
    conf_read_string(conf,c_dum,"nz_filename",par->fnameNzSrcs[ii]);
    conf_read_string(conf,c_dum,"bias_filename",par->fnameBzSrcs[ii]);
    conf_read_bool(conf,c_dum,"include_shear",&(par->shear_srcs[ii]));
    if(par->shear_srcs[ii])
      par->do_lensing=1;
    conf_read_bool(conf,c_dum,"store_skewers",&(par->skw_srcs[ii]));
    if(par->skw_srcs[ii])
      par->do_skewers=1;
  }
  if(par->n_srcs>0)
    par->do_srcs=1;

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
    conf_read_double_array(conf,"kappa","z_out",par->z_kappa_out,&(par->n_kappa),NPLANES_MAX);
    conf_read_int(conf,"kappa","nside",&(par->nside_kappa));
  }

  cset=config_lookup(conf,"isw");
  if(cset!=NULL) {
    par->do_isw=1;
    conf_read_double_array(conf,"isw","z_out",par->z_isw_out,&(par->n_isw),NPLANES_MAX);
    conf_read_int(conf,"isw","nside",&(par->nside_isw));
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

#ifdef _ADD_EXTRA_KAPPA
  if(par->do_kappa) {
    par->need_extra_kappa=my_calloc(par->n_kappa,sizeof(int));
    par->fl_mean_extra_kappa=my_malloc(par->n_kappa*sizeof(flouble *));
    par->cl_extra_kappa=my_malloc(par->n_kappa*sizeof(flouble *));
  }
  if(par->do_isw) {
    par->need_extra_isw=my_calloc(par->n_isw,sizeof(int));
    par->fl_mean_extra_isw=my_malloc(par->n_isw*sizeof(flouble *));
    par->cl_extra_isw=my_malloc(par->n_isw*sizeof(flouble *));
  }
#endif //_ADD_EXTRA_KAPPA

  config_destroy(conf);

  if(par->r2_smooth>0) {
    par->r2_smooth=pow(par->r2_smooth,2);
    par->do_smoothing=1;
  }
  else
    par->do_smoothing=0;

  par->need_beaming=par->do_lensing+par->do_kappa+par->do_isw+par->do_skewers;
  init_fftw(par);

  get_max_memory(par,test_memory+par->just_do_pred);

  cosmo_set(par);
  print_info("\n");

  double dk=2*M_PI/par->l_box;
  print_info("Run parameters: \n");
  print_info("  %.3lf < z < %.3lf\n",par->z_min,par->z_max);
  print_info("  %.3lf < r/(Mpc/h) < %.3lf\n",par->r_min,par->r_max);
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
  if(par->do_skewers)
    print_info("  Some populations will include LoS skewers\n");
  if(par->do_imap)
    print_info("  %d intensity mapping species\n",par->n_imap);
  if(par->do_kappa)
    print_info("  %d lensing source planes\n",par->n_kappa);
  if(par->do_isw)
    print_info("  %d ISW source planes\n",par->n_isw);
  if(par->do_lensing)
    print_info("  Will include lensing shear\n");
  if(!par->need_beaming)
    print_info("  Will NOT need to all-to-all communicate fields\n");
  print_info("\n");

  if(test_memory) {
#ifdef _DEBUG
    fclose(par->f_dbg);
#endif //_DEBUG
    free(par);
    return NULL;
  }

  if(par->just_do_pred)
    return par;

  allocate_fftw(par);

  if(par->do_srcs) {
    par->cats_c=my_malloc(par->n_srcs*sizeof(CatalogCartesian *));
    par->cats=my_malloc(par->n_srcs*sizeof(CatalogCartesian *));
    par->nsources_c_this=my_calloc(par->n_srcs,sizeof(long));
    par->nsources_this=my_calloc(par->n_srcs,sizeof(long));
    for(ii=0;ii<par->n_srcs;ii++) {
      par->cats_c[ii]=NULL;
      par->cats[ii]=NULL;
    }
  }

  if(par->do_imap) {
    par->imap=my_malloc(par->n_imap*sizeof(HealpixShells *));
    for(ii=0;ii<par->n_imap;ii++) {
      FILE *fnu=fopen(par->fnameNuImap[ii],"r");
      if(fnu==NULL) error_open_file(par->fnameNuImap[ii]);
      par->imap[ii]=hp_shell_alloc(par->nside_imap[ii],par->nside_base,linecount(fnu));
      fclose(fnu);
    }
  }

  if(par->do_kappa)
    par->kmap=hp_shell_alloc(par->nside_kappa,par->nside_base,par->n_kappa);

  if(par->do_isw)
    par->pd_map=hp_shell_alloc(par->nside_isw,par->nside_base,par->n_isw);
  
  compute_tracer_cosmo(par);
  
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

void write_imap(ParamCoLoRe *par)
{
  int ipop;
  char fname[256];

  if(NodeThis==0) timer(0);
  print_info("*** Writing intensity maps\n");
  for(ipop=0;ipop<par->n_imap;ipop++) {
    int ir;
    HealpixShells *imap=par->imap[ipop];
    long npx=he_nside2npix(imap->nside);
    flouble *map_write=my_malloc(npx*sizeof(flouble));
    int *map_nadd=my_malloc(npx*sizeof(int));
    print_info(" %d-th species\n",ipop);
    for(ir=0;ir<imap->nr;ir++) {
      long ip;
      long ir_t=ir*imap->num_pix;
      
      //Write local pixels to dummy map
      for(ip=0;ip<npx;ip++) {
	map_write[ip]=0;
	map_nadd[ip]=0;
      }
      sprintf(fname,"!%s_imap_s%d_nu%03d.fits",par->prefixOut,ipop+1,ir);
      for(ip=0;ip<npx;ip++) {
	map_write[ip]+=imap->data[ir_t+ip];
	map_nadd[ip]+=imap->nadd[ir_t+ip];
      }
      
      //Collect all dummy maps
#ifdef _HAVE_MPI
      if(NodeThis==0)
	MPI_Reduce(MPI_IN_PLACE,map_write,npx,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
      else
	MPI_Reduce(map_write   ,NULL     ,npx,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
      if(NodeThis==0)
	MPI_Reduce(MPI_IN_PLACE,map_nadd,npx,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
      else
	MPI_Reduce(map_nadd    ,NULL    ,npx,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
#endif //_HAVE_MPI
      
      for(ip=0;ip<npx;ip++) {
	if(map_nadd[ip]>0)
	  map_write[ip]/=map_nadd[ip];
      }
      
      //Write dummy map
      if(NodeThis==0)
	he_write_healpix_map(&map_write,1,imap->nside,fname,0);
    }
    free(map_write);
    free(map_nadd);
  }
  if(NodeThis==0) timer(2);
  print_info("\n");
}

void write_kappa(ParamCoLoRe *par)
{
  int ir;
  char fname[256];
  long npx=he_nside2npix(par->nside_kappa);
  flouble *map_write=my_malloc(npx*sizeof(flouble));
  int *map_nadd=my_malloc(npx*sizeof(int));
  if(NodeThis==0) timer(0);
  print_info("*** Writing kappa source maps\n");
  for(ir=0;ir<par->kmap->nr;ir++) {
    long ip;
    long ir_t=ir*par->kmap->num_pix;
    
    //Write local pixels to dummy map
    for(ip=0;ip<npx;ip++) {
      map_write[ip]=0;
      map_nadd[ip]=0;
    }
    sprintf(fname,"!%s_kappa_z%03d.fits",par->prefixOut,ir);
    for(ip=0;ip<par->kmap->num_pix;ip++) {
      int id_pix=par->kmap->listpix[ip];
      map_write[id_pix]+=par->kmap->data[ir_t+ip];
      map_nadd[ id_pix]+=par->kmap->nadd[ir_t+ip];
    }
    
    //Collect all dummy maps
#ifdef _HAVE_MPI
    if(NodeThis==0)
      MPI_Reduce(MPI_IN_PLACE,map_write,npx,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    else
      MPI_Reduce(map_write   ,NULL     ,npx,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    if(NodeThis==0)
      MPI_Reduce(MPI_IN_PLACE,map_nadd,npx,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    else
      MPI_Reduce(map_nadd    ,NULL    ,npx,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
#endif //_HAVE_MPI
    
    for(ip=0;ip<npx;ip++) {
      if(map_nadd[ip]>0)
	map_write[ip]/=map_nadd[ip];
    }
    
#ifdef _ADD_EXTRA_KAPPA
    if(par->need_extra_kappa[ir]) {
      if(NodeThis==0) {
	int lmax=3*par->kmap->nside;
	flouble *map_extra;
	flouble *map_mean=my_calloc(npx,sizeof(flouble));
	fcomplex *alm=my_malloc(he_nalms(lmax)*sizeof(fcomplex));
	print_info("Adding perturbations to kappa shell #%d\n",ir+1);
	he_map2alm(par->kmap->nside,lmax,1,&map_write,&alm);
	he_alter_alm(lmax,0,alm,alm,par->fl_mean_extra_kappa[ir]);
	he_alm2map(par->kmap->nside,lmax,1,&map_mean,&alm);
	free(alm);
	map_extra=he_synfast(par->cl_extra_kappa[ir],par->kmap->nside,lmax,par->seed_rng);
	for(ip=0;ip<npx;ip++)
	  map_write[ip]=map_write[ip]+map_mean[ip]+map_extra[ip];
	free(map_extra);
	free(map_mean);
      }
    }
#endif //_ADD_EXTRA_KAPPA
    
    //Write dummy map
    if(NodeThis==0)
      he_write_healpix_map(&map_write,1,par->nside_kappa,fname,1);
  }
  free(map_write);
  free(map_nadd);
  if(NodeThis==0) timer(2);
  print_info("\n");
}

void write_isw(ParamCoLoRe *par)
{
  int ir;
  char fname[256];
  long npx=he_nside2npix(par->nside_isw);
  flouble *map_write=my_malloc(npx*sizeof(flouble));
  int *map_nadd=my_malloc(npx*sizeof(int));
  if(NodeThis==0) timer(0);
  print_info("*** Writing isw source maps\n");
  for(ir=0;ir<par->pd_map->nr;ir++) {
    long ip;
    long ir_t=ir*par->pd_map->num_pix;
    
    //Write local pixels to dummy map
    for(ip=0;ip<npx;ip++) {
      map_write[ip]=0;
      map_nadd[ip]=0;
    }
    sprintf(fname,"!%s_isw_z%03d.fits",par->prefixOut,ir);
    for(ip=0;ip<par->pd_map->num_pix;ip++) {
      int id_pix=par->pd_map->listpix[ip];
      map_write[id_pix]+=par->pd_map->data[ir_t+ip];
      map_nadd[ id_pix]+=par->pd_map->nadd[ir_t+ip];
    }
    
    //Collect all dummy maps
#ifdef _HAVE_MPI
    if(NodeThis==0)
      MPI_Reduce(MPI_IN_PLACE,map_write,npx,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    else
      MPI_Reduce(map_write   ,NULL     ,npx,FLOUBLE_MPI,MPI_SUM,0,MPI_COMM_WORLD);
    if(NodeThis==0)
      MPI_Reduce(MPI_IN_PLACE,map_nadd,npx,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    else
      MPI_Reduce(map_nadd    ,NULL    ,npx,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
#endif //_HAVE_MPI
    
    for(ip=0;ip<npx;ip++) {
      if(map_nadd[ip]>0)
	map_write[ip]/=map_nadd[ip];
    }
    
#ifdef _ADD_EXTRA_KAPPA
    if(par->need_extra_isw[ir]) {
      if(NodeThis==0) {
	int lmax=3*par->pd_map->nside;
	flouble *map_extra;
	flouble *map_mean=my_calloc(npx,sizeof(flouble));
	fcomplex *alm=my_malloc(he_nalms(lmax)*sizeof(fcomplex));
	print_info("Adding perturbations to isw shell #%d\n",ir+1);
	he_map2alm(par->pd_map->nside,lmax,1,&map_write,&alm);
	he_alter_alm(lmax,0,alm,alm,par->fl_mean_extra_isw[ir]);
	he_alm2map(par->pd_map->nside,lmax,1,&map_mean,&alm);
	free(alm);
	map_extra=he_synfast(par->cl_extra_isw[ir],par->pd_map->nside,lmax,par->seed_rng);
	for(ip=0;ip<npx;ip++)
	  map_write[ip]=map_write[ip]+map_mean[ip]+map_extra[ip];
	free(map_extra);
	free(map_mean);
      }
    }
#endif //_ADD_EXTRA_KAPPA
    
    //Write dummy map
    if(NodeThis==0)
      he_write_healpix_map(&map_write,1,par->nside_isw,fname,1);
  }
  free(map_write);
  free(map_nadd);
  if(NodeThis==0) timer(2);
  print_info("\n");
}

static void write_catalog(ParamCoLoRe *par,int ipop)
{
  char fname[256];

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
    
    print_info(" %d-th population (HDF5)\n",ipop);
    sprintf(fname,"%s_srcs_s%d_%d.h5",par->prefixOut,ipop+1,NodeThis);
    
    //Create file
    file_id=H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    //Write table for each galaxy type
    char table_title[256],table_name[256];
    sprintf(table_title,"sources%d_data",ipop+1);
    sprintf(table_name,"/sources%d",ipop+1);
    H5TBmake_table(table_title,file_id,table_name,6,par->nsources_this[ipop],sizeof(Src),
		   names,dst_offset,gal_types,chunk_size,NULL,0,par->cats[ipop]->srcs);
    H5LTset_attribute_string(file_id,table_name,"FIELD_0_UNITS",tunit[0]);
    H5LTset_attribute_string(file_id,table_name,"FIELD_1_UNITS",tunit[1]);
    H5LTset_attribute_string(file_id,table_name,"FIELD_2_UNITS",tunit[2]);
    H5LTset_attribute_string(file_id,table_name,"FIELD_3_UNITS",tunit[3]);
    H5LTset_attribute_string(file_id,table_name,"FIELD_4_UNITS",tunit[4]);
    H5LTset_attribute_string(file_id,table_name,"FIELD_5_UNITS",tunit[5]);
    
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
    
    print_info(" %d-th population (FITS)\n",ipop);
    sprintf(fname,"!%s_srcs_s%d_%d.fits",par->prefixOut,ipop+1,NodeThis);
    
    fits_create_file(&fptr,fname,&status);
    fits_create_tbl(fptr,BINARY_TBL,0,nfields,ttype,tform,tunit,NULL,&status);
    fits_update_key(fptr,TSTRING,"CONTENTS","Source catalog",NULL,&status);
    
    fits_get_rowsize(fptr,&nrw,&status);
    type_arr=my_malloc(nrw*sizeof(int));
    ra_arr=my_malloc(nrw*sizeof(float));
    dec_arr=my_malloc(nrw*sizeof(float));
    z0_arr=my_malloc(nrw*sizeof(float));
    rsd_arr=my_malloc(nrw*sizeof(float));
    e1_arr=my_malloc(nrw*sizeof(float));
    e2_arr=my_malloc(nrw*sizeof(float));
    
    long row_here=0;
    while(row_here<par->nsources_this[ipop]) {
      long nrw_here;
      if(row_here+nrw>par->nsources_this[ipop])
	nrw_here=par->nsources_this[ipop]-row_here;
      else
	nrw_here=nrw;
      
      for(ii=0;ii<nrw_here;ii++) {
	type_arr[ii]=ipop;
	ra_arr[ii]=par->cats[ipop]->srcs[row_here+ii].ra;
	dec_arr[ii]=par->cats[ipop]->srcs[row_here+ii].dec;
	z0_arr[ii]=par->cats[ipop]->srcs[row_here+ii].z0;
	rsd_arr[ii]=par->cats[ipop]->srcs[row_here+ii].dz_rsd;
	if(par->do_lensing) {
	  e1_arr[ii]=par->cats[ipop]->srcs[row_here+ii].e1;
	  e2_arr[ii]=par->cats[ipop]->srcs[row_here+ii].e2;
	}
      }
      fits_write_col(fptr,TINT  ,1,row_here+1,1,nrw_here,type_arr,&status);
      fits_write_col(fptr,TFLOAT,2,row_here+1,1,nrw_here,ra_arr,&status);
      fits_write_col(fptr,TFLOAT,3,row_here+1,1,nrw_here,dec_arr,&status);
      fits_write_col(fptr,TFLOAT,4,row_here+1,1,nrw_here,z0_arr,&status);
      fits_write_col(fptr,TFLOAT,5,row_here+1,1,nrw_here,rsd_arr,&status);
      if(par->do_lensing) {
	fits_write_col(fptr,TFLOAT,6,row_here+1,1,nrw_here,e1_arr,&status);
	fits_write_col(fptr,TFLOAT,7,row_here+1,1,nrw_here,e2_arr,&status);
      }
      
      row_here+=nrw_here;
    }
    
    if(par->cats[ipop]->has_skw) {
      int ir;
      long nelements,naxis=2;
      long naxes[2];
      
      //Write density skewers
      naxes[0]=par->cats[ipop]->nr;
      naxes[1]=par->cats[ipop]->nsrc;
      nelements=naxes[0]*naxes[1];
      fits_create_img(fptr,FLOAT_IMG,naxis,naxes,&status);
      fits_update_key(fptr, TSTRING, "CONTENTS", "density skewers",NULL, &status);
      fits_write_img(fptr,TFLOAT,1,nelements,par->cats[ipop]->d_skw,&status);
      
      //Write velocity skewers
      fits_create_img(fptr,FLOAT_IMG,naxis,naxes,&status);
      fits_update_key(fptr, TSTRING, "CONTENTS", "velocity skewers",NULL, &status);
      fits_write_img(fptr,TFLOAT,1,nelements,par->cats[ipop]->v_skw,&status);
      
      //Write slicing information
      float sg=sqrt(par->sigma2_gauss);
      float *ra=my_malloc(par->cats[ipop]->nr*sizeof(float));
      float *za=my_malloc(par->cats[ipop]->nr*sizeof(float));
      float *gda=my_malloc(par->cats[ipop]->nr*sizeof(float));
      float *gva=my_malloc(par->cats[ipop]->nr*sizeof(float));
      char *tt[]={"R","Z","D","V"};
      char *tf[]={"1E","1E","1E","1E"};
      char *tu[]={"MPC_H","NA","NA","NA"};
      for(ir=0;ir<par->cats[ipop]->nr;ir++) {
	double r=par->cats[ipop]->dr*(ir+0.5);
	ra[ir]=r;
	za[ir]=get_bg(par,r,BG_Z,0);
	gda[ir]=get_bg(par,r,BG_D1,0);
	gva[ir]=get_bg(par,r,BG_V1,0);
      }
      fits_create_tbl(fptr,BINARY_TBL,0,4,tt,tf,tu,NULL,&status);
      fits_update_key(fptr,TSTRING,"CONTENTS","Background cosmology",NULL,&status);
      fits_update_key(fptr,TFLOAT,"SIGMA_G",&sg,NULL, &status);
      fits_write_col(fptr,TFLOAT,1,1,1,par->cats[ipop]->nr,ra,&status);
      fits_write_col(fptr,TFLOAT,2,1,1,par->cats[ipop]->nr,za,&status);
      fits_write_col(fptr,TFLOAT,3,1,1,par->cats[ipop]->nr,gda,&status);
      fits_write_col(fptr,TFLOAT,4,1,1,par->cats[ipop]->nr,gva,&status);
      free(ra); free(za); free(gda); free(gva);
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
  else {
    print_info(" %d-th population (ASCII)\n",ipop);
    sprintf(fname,"%s_srcs_s%d_%d.txt",par->prefixOut,ipop+1,NodeThis);
    
    long jj;
    FILE *fil=fopen(fname,"w");
    if(fil==NULL) error_open_file(fname);
    fprintf(fil,"#[1]type [2]RA, [3]dec, [4]z0, [5]dz_RSD ");
    if(par->do_lensing)
      fprintf(fil,"#[6]e1, [7]e2\n");
    else
      fprintf(fil,"\n");
    for(jj=0;jj<par->nsources_this[ipop];jj++) {
      fprintf(fil,"%d %E %E %E %E ",
	      ipop,par->cats[ipop]->srcs[jj].ra,par->cats[ipop]->srcs[jj].dec,
	      par->cats[ipop]->srcs[jj].z0,par->cats[ipop]->srcs[jj].dz_rsd);
      if(par->do_lensing)
	fprintf(fil,"%E %E \n",par->cats[ipop]->srcs[jj].e1,par->cats[ipop]->srcs[jj].e2);
      else
	fprintf(fil,"\n");
    }
    fclose(fil);
  }
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
    for(ii=0;ii<par->n_srcs;ii++) {
      free(par->srcs_bz_arr[ii]);
      free(par->srcs_nz_arr[ii]);
      free(par->srcs_norm_arr[ii]);
      if(par->cats_c!=NULL) {
	if(par->cats_c[ii]!=NULL)
	  catalog_cartesian_free(par->cats_c[ii]);
      }
      if(par->cats!=NULL) {
	if(par->cats[ii]!=NULL)
	  catalog_free(par->cats[ii]);
      }
    }
    if(par->cats_c!=NULL)
      free(par->cats_c);
    if(par->cats!=NULL)
      free(par->cats);
    free(par->nsources_c_this);
    free(par->nsources_this);
  }

  if(par->do_imap) {
    for(ii=0;ii<par->n_imap;ii++) {
      free(par->imap_bz_arr[ii]);
      free(par->imap_tz_arr[ii]);
      free(par->imap_norm_arr[ii]);
      if(par->imap!=NULL)
	hp_shell_free(par->imap[ii]);
    }
    if(par->imap!=NULL)
      free(par->imap);
  }

  if(par->do_kappa) {
    if(par->kmap!=NULL)
      hp_shell_free(par->kmap);
#ifdef _ADD_EXTRA_KAPPA
    for(ii=0;ii<par->n_kappa;ii++) {
      if(par->need_extra_kappa[ii]) {
	free(par->fl_mean_extra_kappa[ii]);
	free(par->cl_extra_kappa[ii]);
      }
    }
    free(par->need_extra_kappa);
    free(par->fl_mean_extra_kappa);
    free(par->cl_extra_kappa);
#endif //_ADD_EXTRA_KAPPA
  }

  if(par->do_isw) {
    if(par->pd_map!=NULL)
      hp_shell_free(par->pd_map);
#ifdef _ADD_EXTRA_KAPPA
    for(ii=0;ii<par->n_isw;ii++) {
      if(par->need_extra_isw[ii]) {
    	free(par->fl_mean_extra_isw[ii]);
    	free(par->cl_extra_isw[ii]);
      }
    }
    free(par->need_extra_isw);
    free(par->fl_mean_extra_isw);
    free(par->cl_extra_isw);
#endif //_ADD_EXTRA_KAPPA
  }


#ifdef _DEBUG
  fclose(par->f_dbg);
#endif //_DEBUG
}
