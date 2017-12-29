#include "common.h"

//HE_IO
void he_write_healpix_map(flouble **tmap,int nfields,long nside,char *fname,int isnest)
{
  fitsfile *fptr;
  int ii,status=0;
  char **ttype,**tform,**tunit;
  float *map_dum=my_malloc(nside2npix(nside)*sizeof(float));

  ttype=my_malloc(nfields*sizeof(char *));
  tform=my_malloc(nfields*sizeof(char *));
  tunit=my_malloc(nfields*sizeof(char *));
  for(ii=0;ii<nfields;ii++) {
    ttype[ii]=my_malloc(256);
    tform[ii]=my_malloc(256);
    tunit[ii]=my_malloc(256);
    sprintf(ttype[ii],"map %d",ii+1);
    sprintf(tform[ii],"1E");
    sprintf(tunit[ii],"uK");
  }

  fits_create_file(&fptr,fname,&status);
  fits_create_tbl(fptr,BINARY_TBL,0,nfields,ttype,tform,
		  tunit,"BINTABLE",&status);
  fits_write_key(fptr,TSTRING,"PIXTYPE","HEALPIX","HEALPIX Pixelisation",
		 &status);

  fits_write_key(fptr,TSTRING,"ORDERING","RING",
		 "Pixel ordering scheme, either RING or NESTED",&status);
  fits_write_key(fptr,TLONG,"NSIDE",&nside,
		 "Resolution parameter for HEALPIX",&status);
  fits_write_key(fptr,TSTRING,"COORDSYS","G",
		 "Pixelisation coordinate system",&status);
  fits_write_comment(fptr,
		     "G = Galactic, E = ecliptic, C = celestial = equatorial",
		     &status);
  for(ii=0;ii<nfields;ii++) {
    long ip;
    if(isnest)
      he_nest2ring_inplace(tmap[ii],nside);
    for(ip=0;ip<nside2npix(nside);ip++)
      map_dum[ip]=(float)(tmap[ii][ip]);
    fits_write_col(fptr,TFLOAT,ii+1,1,1,nside2npix(nside),map_dum,&status);
  }
  fits_close_file(fptr, &status);

  for(ii=0;ii<nfields;ii++) {
    free(ttype[ii]);
    free(tform[ii]);
    free(tunit[ii]);
  }
  free(ttype);
  free(tform);
  free(tunit);
}

flouble *he_read_healpix_map(char *fname,long *nside,int nfield)
{
  //////
  // Reads a healpix map from file fname. The map will be
  // read from column #nfield. It also returns the map's nside.
  int status=0,hdutype,nfound,anynul;
  long naxes,*naxis,npix;
  fitsfile *fptr;
  flouble *map,nulval;
  char order_in_file[32];
  int nested_in_file=0;

  fits_open_file(&fptr,fname,READONLY,&status);
  fits_movabs_hdu(fptr,2,&hdutype,&status);
  fits_read_key_lng(fptr,"NAXIS",&naxes,NULL,&status);
  naxis=my_malloc(naxes*sizeof(long));
  fits_read_keys_lng(fptr,"NAXIS",1,naxes,naxis,&nfound,&status);
  fits_read_key_lng(fptr,"NSIDE",nside,NULL,&status);
  npix=12*(*nside)*(*nside);
  if(npix%naxis[1]!=0) {
    fprintf(stderr,"WTFFF\n");
    exit(1);
  }

  if (fits_read_key(fptr, TSTRING, "ORDERING", order_in_file, NULL, &status))
    report_error(1,"Could not find %s keyword in in file %s\n","ORDERING",fname);
  if(!strncmp(order_in_file,"NEST",4))
    nested_in_file=1;

  map=my_malloc(npix*sizeof(flouble));
#ifdef _SPREC
  fits_read_col(fptr,TFLOAT,nfield+1,1,1,npix,&nulval,map,&anynul,&status);
#else //_SPREC
  fits_read_col(fptr,TDOUBLE,nfield+1,1,1,npix,&nulval,map,&anynul,&status);
#endif //_SPREC
  free(naxis);

  fits_close_file(fptr,&status);

  flouble *map_ring;
  if(nested_in_file) {
    long ipring,ipnest;

    printf("read_healpix_map: input is nested. Transforming to ring.\n");
    map_ring=my_malloc(npix*sizeof(flouble));
    for(ipnest=0;ipnest<npix;ipnest++) {
      nest2ring(*nside,ipnest,&ipring);
      map_ring[ipring]=map[ipnest];
    }
    free(map);
  }
  else
    map_ring=map;

  return map_ring;
}


//HE_PIX
static double fmodulo (double v1, double v2)
{
  if (v1>=0)
    return (v1<v2) ? v1 : fmod(v1,v2);
  double tmp=fmod(v1,v2)+v2;
  return (tmp==v2) ? 0. : tmp;
}
static int imodulo (int v1, int v2)
{ int v=v1%v2; return (v>=0) ? v : v+v2; }

long he_nside2npix(long nside)
{
  return 12*nside*nside;
}

double he_pixel_area(int nside)
{
  return M_PI/(3*nside*nside);
}

static const double twopi=6.283185307179586476925286766559005768394;
static const double twothird=2.0/3.0;
static const double inv_halfpi=0.6366197723675813430755350534900574;
long he_ang2pix(long nside,double cth,double phi)
{
  double ctha=fabs(cth);
  double tt=fmodulo(phi,twopi)*inv_halfpi; /* in [0,4) */
  
  if (ctha<=twothird) {/* Equatorial region */
    double temp1=nside*(0.5+tt);
    double temp2=nside*cth*0.75;
    int jp=(int)(temp1-temp2); /* index of  ascending edge line */
    int jm=(int)(temp1+temp2); /* index of descending edge line */
    int ir=nside+1+jp-jm; /* ring number counted from cth=2/3 */ /* in {1,2n+1} */
    int kshift=1-(ir&1); /* kshift=1 if ir even, 0 otherwise */
    int ip=(jp+jm-nside+kshift+1)/2; /* in {0,4n-1} */
    ip=imodulo(ip,4*nside);

    return nside*(nside-1)*2 + (ir-1)*4*nside + ip;
  }
  else {  /* North & South polar caps */
    double tp=tt-(int)(tt);
    double tmp=nside*sqrt(3*(1-ctha));
    int jp=(int)(tp*tmp); /* increasing edge line index */
    int jm=(int)((1.0-tp)*tmp); /* decreasing edge line index */
    int ir=jp+jm+1; /* ring number counted from the closest pole */
    int ip=(int)(tt*ir); /* in {0,4*ir-1} */
    ip = imodulo(ip,4*ir);

    if (cth>0)
      return 2*ir*(ir-1)+ip;
    else
      return 12*nside*nside-2*ir*(ir+1)+ip;
  }
}

int he_ring_num(long nside,double z)
{
  //Returns ring index for normalized height z
  int iring;

  iring=(int)(nside*(2-1.5*z)+0.5);
  if(z>0.66666666) {
    iring=(int)(nside*sqrt(3*(1-z))+0.5);
    if(iring==0) iring=1;
  }

  if(z<-0.66666666) {
    iring=(int)(nside*sqrt(3*(1+z))+0.5);
    if(iring==0) iring=1;
    iring=4*nside-iring;
  }

  return iring;
}

static void get_ring_limits(long nside,int iz,long *ip_lo,long *ip_hi)
{
  long ir;
  long ipix1,ipix2;
  long npix=12*nside*nside;
  long ncap=2*nside*(nside-1);

  if((iz>=nside)&&(iz<=3*nside)) { //eqt
    ir=iz-nside+1;
    ipix1=ncap+4*nside*(ir-1);
    ipix2=ipix1+4*nside-1;
  }
  else {
    if(iz<nside) { //north
      ir=iz;
      ipix1=2*ir*(ir-1);
      ipix2=ipix1+4*ir-1;
    }
    else { //south
      ir=4*nside-iz;
      ipix1=npix-2*ir*(ir+1);
      ipix2=ipix1+4*ir-1;
    }
  }

  *ip_lo=ipix1;
  *ip_hi=ipix2;
}

long *he_query_strip(long nside,double theta1,double theta2,
		     long *npix_strip)
{
  long *pixlist;
  double z_hi=cos(theta1);
  double z_lo=cos(theta2);
  int irmin,irmax;

  if((theta2<=theta1)||
     (theta1<0)||(theta1>M_PI)||
     (theta2<0)||(theta2>M_PI)) {
    report_error(1,"wrong strip boundaries\n");
  }

  irmin=he_ring_num(nside,z_hi);
  irmax=he_ring_num(nside,z_lo);

  //Count number of pixels in strip
  int iz;
  long npix_in_strip=0;
  for(iz=irmin;iz<=irmax;iz++) {
    long ipix1,ipix2;
    get_ring_limits(nside,iz,&ipix1,&ipix2);
    npix_in_strip+=ipix2-ipix1+1;
  }
  *npix_strip=npix_in_strip;
  pixlist=my_malloc(npix_in_strip*sizeof(long));

  //Count number of pixels in strip
  long i_list=0;
  for(iz=irmin;iz<=irmax;iz++) {
    long ipix1,ipix2,ip;
    get_ring_limits(nside,iz,&ipix1,&ipix2);
    for(ip=ipix1;ip<=ipix2;ip++) {
      pixlist[i_list]=ip;
      i_list++;
    }    
  }

  return pixlist;
}

void he_ring2nest_inplace(flouble *map_in,long nside)
{
  long npix=12*nside*nside;
  flouble *map_out=my_malloc(npix*sizeof(flouble));

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(map_in,nside,npix,map_out)
#endif //_HAVE_OMP
  {
    long ip;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ip=0;ip<npix;ip++) {
      long inest;
      ring2nest(nside,ip,&inest);
      
      map_out[inest]=map_in[ip];
    } //end omp for
  } //end omp parallel
  memcpy(map_in,map_out,npix*sizeof(flouble));
  
  free(map_out);
}

void he_nest2ring_inplace(flouble *map_in,long nside)
{
  long npix=12*nside*nside;
  flouble *map_out=my_malloc(npix*sizeof(flouble));

#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(map_in,nside,npix,map_out)
#endif //_HAVE_OMP
  {
    long ip;

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ip=0;ip<npix;ip++) {
      long iring;
      nest2ring(nside,ip,&iring);

      map_out[iring]=map_in[ip];
    } //end omp for
  } //end omp parallel
  memcpy(map_in,map_out,npix*sizeof(flouble));

  free(map_out);
}

void he_udgrade(flouble *map_in,long nside_in,
		flouble *map_out,long nside_out,
		int nest)
{
  long npix_in=nside2npix(nside_in);
  long npix_out=nside2npix(nside_out);

  if(nside_in>nside_out) {
    long ii;
    long np_ratio=npix_in/npix_out;
    double i_np_ratio=1./((double)np_ratio);
    
    for(ii=0;ii<npix_out;ii++) {
      int jj;
      double tot=0;

      if(nest) {
	for(jj=0;jj<np_ratio;jj++)
	  tot+=map_in[jj+ii*np_ratio];
	map_out[ii]=tot*i_np_ratio;
      }
      else {
	long inest_out;

	ring2nest(nside_out,ii,&inest_out);
	for(jj=0;jj<np_ratio;jj++) {
	  long iring_in;
	  
	  nest2ring(nside_in,jj+np_ratio*inest_out,&iring_in);
	  tot+=map_in[iring_in];
	}
	map_out[ii]=tot*i_np_ratio;
      }
    }
  }
  else {
    long ii;
    long np_ratio=npix_out/npix_in;
    
    for(ii=0;ii<npix_in;ii++) {
      int jj;
      
      if(nest) {
	flouble value=map_in[ii];

	for(jj=0;jj<np_ratio;jj++)
	  map_out[jj+ii*np_ratio]=value;
      }
      else {
	long inest_in;
	flouble value=map_in[ii];
	ring2nest(nside_in,ii,&inest_in);
	
	for(jj=0;jj<np_ratio;jj++) {
	  long iring_out;
	  
	  nest2ring(nside_out,jj+inest_in*np_ratio,&iring_out);
	  map_out[iring_out]=value;
	}
      }
    }
  }
}


#ifdef _WITH_SHT
//HE_SHT
long he_nalms(int lmax)
{
  return ((lmax+1)*(lmax+2))/2;
}

long he_indexlm(int l,int m,int lmax)
{
  if(m>0)
    return l+m*lmax-(m*(m-1))/2;
  else
    return l;
}

static void sht_wrapper(int spin,int lmax,int nside,int ntrans,flouble **maps,fcomplex **alms,int alm2map)
{
  double time=0;
  sharp_alm_info *alm_info;
  sharp_geom_info *geom_info;

  sharp_make_triangular_alm_info(lmax,lmax,1,&alm_info);
  sharp_make_weighted_healpix_geom_info(nside,1,NULL,&geom_info);
  sharp_execute(alm2map,spin,alms,maps,geom_info,alm_info,ntrans,SHT_TYPE,&time,NULL);
  sharp_destroy_geom_info(geom_info);
  sharp_destroy_alm_info(alm_info);
  //  printf("  Took %lf s according to libsharp\n",time);
}

void he_alm2map(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms)
{
  int nbatches,nodd,itrans;
  nbatches=ntrans/HE_MAX_SHT;
  nodd=ntrans%HE_MAX_SHT;

  for(itrans=0;itrans<nbatches;itrans++)
    sht_wrapper(0,lmax,nside,HE_MAX_SHT,&(maps[itrans*HE_MAX_SHT]),&(alms[itrans*HE_MAX_SHT]),SHARP_ALM2MAP);
  if(nodd>0)
    sht_wrapper(0,lmax,nside,nodd,&(maps[nbatches*HE_MAX_SHT]),&(alms[nbatches*HE_MAX_SHT]),SHARP_ALM2MAP);
}

void he_map2alm(int nside,int lmax,int ntrans,flouble **maps,fcomplex **alms)
{
  int nbatches,nodd,itrans;
  nbatches=ntrans/HE_MAX_SHT;
  nodd=ntrans%HE_MAX_SHT;

  for(itrans=0;itrans<nbatches;itrans++)
    sht_wrapper(0,lmax,nside,HE_MAX_SHT,&(maps[itrans*HE_MAX_SHT]),&(alms[itrans*HE_MAX_SHT]),SHARP_MAP2ALM);
  if(nodd>0)
    sht_wrapper(0,lmax,nside,nodd,&(maps[nbatches*HE_MAX_SHT]),&(alms[nbatches*HE_MAX_SHT]),SHARP_MAP2ALM);
}

void he_alm2cl(fcomplex **alms_1,fcomplex **alms_2,
	       int nmaps_1,int nmaps_2,
	       int pol_1,int pol_2,
	       flouble **cls,int lmax)
{
  int i1,index_cl;

  index_cl=0;
  for(i1=0;i1<nmaps_1;i1++) {
    int i2;
    fcomplex *alm1=alms_1[i1];
    for(i2=0;i2<nmaps_2;i2++) {
      int l;
      fcomplex *alm2=alms_2[i2];
      for(l=0;l<=lmax;l++) {
	int m;
	cls[index_cl][l]=creal(alm1[he_indexlm(l,0,lmax)])*creal(alm2[he_indexlm(l,0,lmax)]);

	for(m=1;m<=l;m++) {
	  long index_lm=he_indexlm(l,m,lmax);
	  cls[index_cl][l]+=2*(creal(alm1[index_lm])*creal(alm2[index_lm])+
			       cimag(alm1[index_lm])*cimag(alm2[index_lm]));
	}
	cls[index_cl][l]/=(2*l+1.);
      }
      index_cl++;
    }
  }
}

void he_anafast(flouble **maps_1,flouble **maps_2,
		int nmaps_1,int nmaps_2,
		int pol_1,int pol_2,
		flouble **cls,int nside,int lmax)
{
  fcomplex **alms_1,**alms_2;
  int i1,spin,lmax_here=3*nside-1;

  alms_1=my_malloc(nmaps_1*sizeof(fcomplex *));
  for(i1=0;i1<nmaps_1;i1++)
    alms_1[i1]=my_malloc(he_nalms(lmax_here)*sizeof(fcomplex));
  if(pol_1) {
    if(nmaps_1!=2)
      report_error(1,"Must provide 2 maps for polarization\n");
    spin=2;
  }
  else 
    spin=0;
  sht_wrapper(spin,lmax,nside,1,maps_1,alms_1,SHARP_MAP2ALM);

  if(maps_1==maps_2)
    alms_2=alms_1;
  else {
    alms_2=my_malloc(nmaps_2*sizeof(fcomplex *));
    for(i1=0;i1<nmaps_2;i1++)
      alms_2[i1]=my_malloc(he_nalms(lmax_here)*sizeof(fcomplex));
    if(pol_2) {
      if(nmaps_2!=2)
	report_error(1,"Must provide 2 maps for polarization\n");
      
      spin=2;
    }
    else 
      spin=0;
    sht_wrapper(spin,lmax,nside,1,maps_2,alms_2,SHARP_MAP2ALM);
  }

  he_alm2cl(alms_1,alms_2,nmaps_1,nmaps_2,pol_1,pol_2,cls,lmax);

  for(i1=0;i1<nmaps_1;i1++)
    free(alms_1[i1]);
  free(alms_1);
  if(alms_1!=alms_2) {
    for(i1=0;i1<nmaps_2;i1++)
      free(alms_2[i1]);
    free(alms_2);
  }
}

flouble *he_generate_beam_window(int lmax,flouble fwhm_amin)
{
  long l;
  flouble sigma=HE_FWHM2SIGMA*fwhm_amin;
  flouble *beam=my_malloc((lmax+1)*sizeof(flouble));

  for(l=0;l<=lmax;l++)
    beam[l]=exp(-0.5*l*(l+1)*sigma*sigma);
  
  return beam;
}

void he_alter_alm(int lmax,flouble fwhm_amin,fcomplex *alm_in,
		  fcomplex *alm_out,flouble *window)
{
  flouble *beam;
  int mm;

  if(window==NULL) beam=he_generate_beam_window(lmax,fwhm_amin);
  else beam=window;

  for(mm=0;mm<=lmax;mm++) {
    int ll;
    for(ll=mm;ll<=lmax;ll++) {
      long index=he_indexlm(ll,mm,lmax);

      alm_out[index]=alm_in[index]*beam[ll];
    }
  }

  if(window==NULL)
    free(beam);
}

flouble *he_synfast(flouble *cl,int nside,int lmax,unsigned int seed)
{
  fcomplex *alms;
  int lmax_here=lmax;
  long npix=12*((long)nside)*nside;
  flouble *map=my_malloc(npix*sizeof(flouble));
  gsl_rng *rng=init_rng(seed);

  if(lmax>3*nside-1)
    lmax_here=3*nside-1;
  alms=my_malloc(he_nalms(lmax_here)*sizeof(fcomplex));

  int ll;
  for(ll=0;ll<=lmax_here;ll++) {
    int mm;
    flouble sigma=sqrt(0.5*cl[ll]);
    double r1,r2;
    rng_gauss(rng,&r1,&r2);
    alms[he_indexlm(ll,0,lmax_here)]=(fcomplex)(M_SQRT2*sigma*r1);

    for(mm=1;mm<=ll;mm++) {
      rng_gauss(rng,&r1,&r2);
      alms[he_indexlm(ll,mm,lmax_here)]=(fcomplex)(sigma*(r1+I*r2));
    }
  }
  he_alm2map(nside,lmax_here,1,&map,&alms);
  free(alms);
  end_rng(rng);

  return map;
}


#ifdef _WITH_NEEDLET
//HE_NT
static double func_fx(double x,void *pars)
{
  return exp(-1./(1.-x*x));
}

static double func_psix(double x,gsl_integration_workspace *w,const gsl_function *f)
{
  double result,eresult;
  gsl_integration_qag(f,-1,x,0,HE_NL_INTPREC,1000,GSL_INTEG_GAUSS61,
		      w,&result,&eresult);

  return result*HE_NORM_FT;
}

static double func_phix(double x,double invB,
			gsl_integration_workspace *w,const gsl_function *f)
{
  if(x<0)
    report_error(1,"Something went wrong");
  else if(x<=invB)
    return 1.;
  else if(x<=1)
    return func_psix(1-2.*(x-invB)/(1-invB),w,f);
  else
    return 0.;

  return -1.;
}

void he_nt_end(HE_nt_param *par)
{
  int j;

  gsl_spline_free(par->b_spline);
  gsl_interp_accel_free(par->b_intacc);
  free(par->nside_arr);
  free(par->lmax_arr);
  for(j=0;j<par->nj;j++)
    free(par->b_arr[j]);
  free(par->b_arr);
  free(par);
}

HE_nt_param *he_nt_init(flouble b_nt,int nside0)
{
  HE_nt_param *par=my_malloc(sizeof(HE_nt_param));
  par->b=b_nt;
  par->inv_b=1./b_nt;
  par->b_spline=gsl_spline_alloc(gsl_interp_cspline,HE_NBAND_NX);
  par->b_intacc=gsl_interp_accel_alloc();
  par->nside0=nside0;

  int ii;
  gsl_function F;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  double *xarr=my_malloc(HE_NBAND_NX*sizeof(double));
  double *barr=my_malloc(HE_NBAND_NX*sizeof(double));
  F.function=&func_fx;
  F.params=NULL;
  for(ii=0;ii<HE_NBAND_NX;ii++) {
    xarr[ii]=pow(b_nt,-1.+2.*(ii+0.)/(HE_NBAND_NX-1));
    barr[ii]=sqrt(func_phix(xarr[ii]/b_nt,1./b_nt,w,&F)-
		  func_phix(xarr[ii],1./b_nt,w,&F));
  }
  gsl_spline_init(par->b_spline,xarr,barr,HE_NBAND_NX);
  gsl_integration_workspace_free(w);
  free(xarr);
  free(barr);

  double lmax=3*nside0-1;
  double djmax=log(lmax)/log(b_nt);
  int jmax=(int)(djmax)+1;
  par->nj=jmax+1;
  par->nside_arr=my_malloc(par->nj*sizeof(int));
  par->lmax_arr=my_malloc(par->nj*sizeof(int));
  for(ii=0;ii<par->nj;ii++) {
    double dlmx=pow(par->b,ii+1);
    //    double dlmn=pow(par->b,ii-1);
    int lmx=(int)(dlmx)+1;
    int ns=pow(2,(int)(log((double)lmx)/log(2.))+1);
    par->nside_arr[ii]=MAX((MIN(ns,par->nside0)),HE_NT_NSIDE_MIN);
    par->lmax_arr[ii]=3*par->nside_arr[ii]-1;
  }

  par->b_arr=my_malloc(par->nj*sizeof(flouble *));
  for(ii=0;ii<par->nj;ii++) {
    par->b_arr[ii]=my_calloc(3*nside0,sizeof(flouble));
    he_nt_get_window(par,ii,par->b_arr[ii]);
  }

  int lmx0=(int)(par->b);
  for(ii=0;ii<=lmx0;ii++) {
    flouble b1=par->b_arr[0][ii];
    flouble b2=par->b_arr[1][ii];
    flouble h_here=b1*b1+b2*b2;
    if(h_here<1)
      par->b_arr[0][ii]=sqrt(1-b2*b2);
  }

  return par;
}

flouble ***he_alloc_needlet(HE_nt_param *par,int pol)
{
  int ii,nmaps;
  flouble ***nt=my_malloc(par->nj*sizeof(flouble **));
  if(pol)
    nmaps=3;
  else
    nmaps=1;

  for(ii=0;ii<par->nj;ii++) {
    int imap;
    long ns=par->nside_arr[ii];

    nt[ii]=my_malloc(nmaps*sizeof(flouble *));
    for(imap=0;imap<nmaps;imap++)
      nt[ii][imap]=my_malloc(12*ns*ns*sizeof(flouble));
  }

  return nt;
}

void he_free_needlet(HE_nt_param *par,int pol,flouble ***nt)
{
  int ii,nmaps;

  if(pol)
    nmaps=3;
  else
    nmaps=1;

  for(ii=0;ii<par->nj;ii++) {
    int imap;
    for(imap=0;imap<nmaps;imap++) {
      free(nt[ii][imap]);
    }
    free(nt[ii]);
  }
  free(nt);
}

void he_nt_get_window(HE_nt_param *par,int j,flouble *b)
{
  int l;
  int lmx=par->lmax_arr[j];
  flouble bfac=1./pow(par->b,j);
  
  for(l=0;l<=lmx;l++) {
    double x=(double)l*bfac;
    if((x<par->inv_b)||(x>par->b))
      b[l]=0;
    else
      b[l]=gsl_spline_eval(par->b_spline,x,par->b_intacc);
  }
}

fcomplex **he_needlet2map(HE_nt_param *par,flouble **map,flouble ***nt,
			  int return_alm,int pol,int qu_in,int qu_out)
{
  int j,nmaps;
  fcomplex **alm,**alm_dum;
  int l_max=3*par->nside0-1;
  long n_alms=he_nalms(l_max);

  //Figure out spin
  if(pol)
    nmaps=3;
  else
    nmaps=1;

  //Allocate alms
  alm=my_malloc(nmaps*sizeof(fcomplex *));
  alm_dum=my_malloc(nmaps*sizeof(fcomplex *));
  for(j=0;j<nmaps;j++) {
    alm[j]=my_calloc(n_alms,sizeof(fcomplex));
    alm_dum[j]=my_malloc(n_alms*sizeof(fcomplex));
  }

  //Loop over scales
  for(j=0;j<par->nj;j++) {
    int mm;
    int lmx=par->lmax_arr[j];
    int imap;

    //Compute alm's for j-th needlet
    sht_wrapper(0,par->lmax_arr[j],par->nside_arr[j],1,nt[j],alm_dum,SHARP_MAP2ALM);
    if(pol) {
      if(qu_in)
	sht_wrapper(2,par->lmax_arr[j],par->nside_arr[j],1,&(nt[j][1]),&(alm_dum[1]),SHARP_MAP2ALM);
      else
	sht_wrapper(0,par->lmax_arr[j],par->nside_arr[j],2,&(nt[j][1]),&(alm_dum[1]),SHARP_MAP2ALM);
    }
    //Loop over spin components
    for(imap=0;imap<nmaps;imap++) {
      //Multiply by window and add to total alm
      for(mm=0;mm<=lmx;mm++) {
	int ll;
	for(ll=mm;ll<=lmx;ll++) {
	  long index0=he_indexlm(ll,mm,l_max);
	  long index=he_indexlm(ll,mm,lmx);
	  alm[imap][index0]+=par->b_arr[j][ll]*alm_dum[imap][index];
	}
      }
    }
  }

  //Transform total alm back to map
  sht_wrapper(0,l_max,par->nside0,1,map,alm,SHARP_ALM2MAP);
  if(pol) {
    if(qu_out)
      sht_wrapper(2,l_max,par->nside0,1,&(map[1]),&(alm[1]),SHARP_ALM2MAP);
    else
      sht_wrapper(0,l_max,par->nside0,2,&(map[1]),&(alm[1]),SHARP_ALM2MAP);
  }

  if(!return_alm) {
    for(j=0;j<nmaps;j++)
      free(alm[j]);
    free(alm);
    alm=NULL;
  }
  for(j=0;j<nmaps;j++)
    free(alm_dum[j]);
  free(alm_dum);

  return alm;
}

fcomplex **he_map2needlet(HE_nt_param *par,flouble **map,flouble ***nt,
			  int return_alm,int pol,int qu_in,int qu_out)
{
  int j,nmaps;
  fcomplex **alm,**alm_dum;
  int l_max=3*par->nside0-1;
  long n_alms=he_nalms(l_max);

  //Figure out spin
  if(pol)
    nmaps=3;
  else
    nmaps=1;

  //Allocate alms
  alm=my_malloc(nmaps*sizeof(fcomplex *));
  alm_dum=my_malloc(nmaps*sizeof(fcomplex *));
  for(j=0;j<nmaps;j++) {
    alm[j]=my_malloc(n_alms*sizeof(fcomplex));
    alm_dum[j]=my_malloc(n_alms*sizeof(fcomplex));
  }

  //SHT
  sht_wrapper(0,l_max,par->nside0,1,map,alm,SHARP_MAP2ALM);
  if(qu_in)
    sht_wrapper(2,l_max,par->nside0,1,&(map[1]),&(alm[1]),SHARP_MAP2ALM);
  else
    sht_wrapper(0,l_max,par->nside0,2,&(map[1]),&(alm[1]),SHARP_MAP2ALM);

  //Iterate over scales
  for(j=0;j<par->nj;j++) {
    int mm,imap;
    int lmx=par->lmax_arr[j];

    //Loop over spin components
    for(imap=0;imap<nmaps;imap++) {
      //Set alms and window to zero
      memset(alm_dum[imap],0,n_alms*sizeof(fcomplex));

      //Multiply alms by window function
      for(mm=0;mm<=lmx;mm++) {
	int ll;
	for(ll=mm;ll<=lmx;ll++) {
	  long index0=he_indexlm(ll,mm,l_max);
	  long index=he_indexlm(ll,mm,lmx);
	  alm_dum[imap][index]=par->b_arr[j][ll]*alm[imap][index0];
	}
      }
    }

    //SHT^-1
    sht_wrapper(0,par->lmax_arr[j],par->nside_arr[j],1,nt[j],alm_dum,SHARP_ALM2MAP);
    if(qu_out)
      sht_wrapper(2,par->lmax_arr[j],par->nside_arr[j],1,&(nt[j][1]),&(alm_dum[1]),SHARP_ALM2MAP);
    else
      sht_wrapper(0,par->lmax_arr[j],par->nside_arr[j],2,&(nt[j][1]),&(alm_dum[1]),SHARP_ALM2MAP);
  }

  if(!return_alm) {
    for(j=0;j<nmaps;j++)
      free(alm[j]);
    free(alm);
    alm=NULL;
  }
  for(j=0;j<nmaps;j++)
    free(alm_dum[j]);
  free(alm_dum);
  
  return alm;
}
#endif //_WITH_NEEDLET
#endif //_WITH_SHT
