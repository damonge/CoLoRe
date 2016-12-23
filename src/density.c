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

static dftw_complex *dftw_alloc_complex(ptrdiff_t dsize)
{
  dftw_complex *p;
#ifdef _SPREC
  p=fftwf_alloc_complex(dsize);
#else //_SPREC
  p=fftw_alloc_complex(dsize);
#endif //_SPREC
  if(p==NULL)
    report_error(1,"Ran out of memory\n");
  return p;
}

static void pos_2_ngp(ParamCoLoRe *par,unsigned long long np,
		      flouble *x,flouble *y,flouble *z,flouble *delta)
{
  unsigned long long ii;
  flouble i_agrid=par->n_grid/par->l_box;
  lint ngx=2*(par->n_grid/2+1);

  for(ii=0;ii<np;ii++) {
    lint index;
    int ax,i0[3];

    i0[0]=(int)(x[ii]*i_agrid+0.5);
    i0[1]=(int)(y[ii]*i_agrid+0.5);
    i0[2]=(int)(z[ii]*i_agrid+0.5);
    for(ax=0;ax<3;ax++) {
      if(i0[ax]>=par->n_grid) i0[ax]-=par->n_grid;
      if(i0[ax]<0) i0[ax]+=par->n_grid;
    }
    i0[2]-=par->iz0_here;

    if((i0[2]>=0) && (i0[2]<par->nz_here)) {
      index=i0[0]+par->n_grid*(i0[1]+ngx*i0[0]);
      delta[index]+=1.;
    }
  }
}

static void pos_2_cic(ParamCoLoRe *par,unsigned long long np,
		      flouble *x,flouble *y,flouble *z,flouble *delta)
{
  unsigned long long ii;
  flouble i_agrid=par->n_grid/par->l_box;
  lint ngx=2*(par->n_grid/2+1);

  for(ii=0;ii<np;ii++) {
    int ax,i0[3],i1[3];
    flouble a0[3],a1[3];

    i0[0]=(int)(x[ii]*i_agrid);
    i0[1]=(int)(y[ii]*i_agrid);
    i0[2]=(int)(z[ii]*i_agrid);
    a1[0]=x[ii]*i_agrid-i0[0];
    a1[1]=y[ii]*i_agrid-i0[1];
    a1[2]=z[ii]*i_agrid-i0[2];
    for(ax=0;ax<3;ax++) {
      a0[ax]=1-a1[ax];
      i1[ax]=i0[ax]+1;
      if(i0[ax]<0) i0[ax]+=par->n_grid;
      if(i1[ax]<0) i1[ax]+=par->n_grid;
      if(i0[ax]>=par->n_grid) i0[ax]-=par->n_grid;
      if(i1[ax]>=par->n_grid) i1[ax]-=par->n_grid;
    }
    i0[2]-=par->iz0_here;
    i1[2]-=par->iz0_here;

    if((i0[2]>=0) && (i0[2]<par->nz_here)) {
      delta[i0[0]+par->n_grid*(i0[1]+ngx*i0[2])]+=a0[0]*a0[1]*a0[2];
      delta[i1[0]+par->n_grid*(i0[1]+ngx*i0[2])]+=a1[0]*a0[1]*a0[2];
      delta[i0[0]+par->n_grid*(i1[1]+ngx*i0[2])]+=a0[0]*a1[1]*a0[2];
      delta[i1[0]+par->n_grid*(i1[1]+ngx*i0[2])]+=a1[0]*a1[1]*a0[2];
    }
    if((i1[2]>=0) && (i1[2]<par->nz_here)) {
      delta[i0[0]+par->n_grid*(i0[1]+ngx*i1[2])]+=a0[0]*a0[1]*a1[2];
      delta[i1[0]+par->n_grid*(i0[1]+ngx*i1[2])]+=a1[0]*a0[1]*a1[2];
      delta[i0[0]+par->n_grid*(i1[1]+ngx*i1[2])]+=a0[0]*a1[1]*a1[2];
      delta[i1[0]+par->n_grid*(i1[1]+ngx*i1[2])]+=a1[0]*a1[1]*a1[2];
    }
  }
}

static void pos_2_tsc(ParamCoLoRe *par,unsigned long long np,
		      flouble *x,flouble *y,flouble *z,flouble *delta)
{
  unsigned long long ii;
  flouble i_agrid=par->n_grid/par->l_box;
  lint ngx=2*(par->n_grid/2+1);

  for(ii=0;ii<np;ii++) {
    int ax,i0[3],ip[3],im[3];
    flouble a0[3],ap[3],am[3];

    i0[0]=(int)(x[ii]*i_agrid+0.5);
    i0[1]=(int)(y[ii]*i_agrid+0.5);
    i0[2]=(int)(z[ii]*i_agrid+0.5);
    a0[0]=x[ii]*i_agrid-i0[0];
    a0[1]=y[ii]*i_agrid-i0[1];
    a0[2]=z[ii]*i_agrid-i0[2];
    for(ax=0;ax<3;ax++) {
      am[ax]=0.5*(0.5-a0[ax])*(0.5-a0[ax]);
      ap[ax]=0.5*(0.5+a0[ax])*(0.5+a0[ax]);
      a0[ax]=0.75-a0[ax]*a0[ax];
      ip[ax]=i0[ax]+1;
      im[ax]=i0[ax]-1;
      if(im[ax]<0) im[ax]+=par->n_grid;
      if(i0[ax]<0) i0[ax]+=par->n_grid;
      if(ip[ax]<0) ip[ax]+=par->n_grid;
      if(im[ax]>=par->n_grid) im[ax]-=par->n_grid;
      if(i0[ax]>=par->n_grid) i0[ax]-=par->n_grid;
      if(ip[ax]>=par->n_grid) ip[ax]-=par->n_grid;
    }
    i0[2]-=par->iz0_here;
    im[2]-=par->iz0_here;
    ip[2]-=par->iz0_here;

    if((im[2]>=0) && (im[2]<par->nz_here)) {
      delta[im[0]+par->n_grid*(im[1]+ngx*im[2])]+=am[0]*am[1]*am[2];
      delta[i0[0]+par->n_grid*(im[1]+ngx*im[2])]+=a0[0]*am[1]*am[2];
      delta[ip[0]+par->n_grid*(im[1]+ngx*im[2])]+=ap[0]*am[1]*am[2];
      delta[im[0]+par->n_grid*(i0[1]+ngx*im[2])]+=am[0]*a0[1]*am[2];
      delta[i0[0]+par->n_grid*(i0[1]+ngx*im[2])]+=a0[0]*a0[1]*am[2];
      delta[ip[0]+par->n_grid*(i0[1]+ngx*im[2])]+=ap[0]*a0[1]*am[2];
      delta[im[0]+par->n_grid*(ip[1]+ngx*im[2])]+=am[0]*ap[1]*am[2];
      delta[i0[0]+par->n_grid*(ip[1]+ngx*im[2])]+=a0[0]*ap[1]*am[2];
      delta[ip[0]+par->n_grid*(ip[1]+ngx*im[2])]+=ap[0]*ap[1]*am[2];
    }
    if((i0[2]>=0) && (i0[2]<par->nz_here)) {
      delta[im[0]+par->n_grid*(im[1]+ngx*i0[2])]+=am[0]*am[1]*a0[2];
      delta[i0[0]+par->n_grid*(im[1]+ngx*i0[2])]+=a0[0]*am[1]*a0[2];
      delta[ip[0]+par->n_grid*(im[1]+ngx*i0[2])]+=ap[0]*am[1]*a0[2];
      delta[im[0]+par->n_grid*(i0[1]+ngx*i0[2])]+=am[0]*a0[1]*a0[2];
      delta[i0[0]+par->n_grid*(i0[1]+ngx*i0[2])]+=a0[0]*a0[1]*a0[2];
      delta[ip[0]+par->n_grid*(i0[1]+ngx*i0[2])]+=ap[0]*a0[1]*a0[2];
      delta[im[0]+par->n_grid*(ip[1]+ngx*i0[2])]+=am[0]*ap[1]*a0[2];
      delta[i0[0]+par->n_grid*(ip[1]+ngx*i0[2])]+=a0[0]*ap[1]*a0[2];
      delta[ip[0]+par->n_grid*(ip[1]+ngx*i0[2])]+=ap[0]*ap[1]*a0[2];
    }
    if((ip[2]>=0) && (ip[2]<par->nz_here)) {
      delta[im[0]+par->n_grid*(im[1]+ngx*ip[2])]+=am[0]*am[1]*ap[2];
      delta[i0[0]+par->n_grid*(im[1]+ngx*ip[2])]+=a0[0]*am[1]*ap[2];
      delta[ip[0]+par->n_grid*(im[1]+ngx*ip[2])]+=ap[0]*am[1]*ap[2];
      delta[im[0]+par->n_grid*(i0[1]+ngx*ip[2])]+=am[0]*a0[1]*ap[2];
      delta[i0[0]+par->n_grid*(i0[1]+ngx*ip[2])]+=a0[0]*a0[1]*ap[2];
      delta[ip[0]+par->n_grid*(i0[1]+ngx*ip[2])]+=ap[0]*a0[1]*ap[2];
      delta[im[0]+par->n_grid*(ip[1]+ngx*ip[2])]+=am[0]*ap[1]*ap[2];
      delta[i0[0]+par->n_grid*(ip[1]+ngx*ip[2])]+=a0[0]*ap[1]*ap[2];
      delta[ip[0]+par->n_grid*(ip[1]+ngx*ip[2])]+=ap[0]*ap[1]*ap[2];
    }
  }
}

static void pos_2_dens(ParamCoLoRe *par,unsigned long long np_here,
		       flouble *x,flouble *y,flouble *z,flouble *dens)
{
  if(par->lpt_interp_type==INTERP_NGP)
    pos_2_ngp(par,np_here,x,y,z,dens);
  else if(par->lpt_interp_type==INTERP_CIC)
    pos_2_cic(par,np_here,x,y,z,dens);
  else if(par->lpt_interp_type==INTERP_TSC)
    pos_2_tsc(par,np_here,x,y,z,dens);
  else
    report_error(1,"Wrong interpolation type\n");
}

static void share_particles(ParamCoLoRe *par,unsigned long long np_allocated,unsigned long long np_real,
			    flouble *x,flouble *y,flouble *z,unsigned long long *np_here)
{
  int inode;
  unsigned long long ip;
  flouble dx=par->l_box/par->n_grid;
  flouble *z_left=my_malloc(NNodes*sizeof(flouble));
  flouble *z_right=my_malloc(NNodes*sizeof(flouble));
  flouble *z_true_left=my_malloc(NNodes*sizeof(flouble));
  flouble *z_true_right=my_malloc(NNodes*sizeof(flouble));
  flouble *z_bleft_left=my_malloc(NNodes*sizeof(flouble));
  flouble *z_bleft_right=my_malloc(NNodes*sizeof(flouble));
  flouble *z_bright_left=my_malloc(NNodes*sizeof(flouble));
  flouble *z_bright_right=my_malloc(NNodes*sizeof(flouble));
  unsigned long long nbuffer=(unsigned long long)(par->lpt_buffer_fraction*par->nz_here*
						  ((lint)(par->n_grid*par->n_grid)));
  flouble *x_b=my_malloc(nbuffer*sizeof(flouble));
  flouble *y_b=my_malloc(nbuffer*sizeof(flouble));
  flouble *z_b=my_malloc(nbuffer*sizeof(flouble));

  for(inode=0;inode<NNodes;inode++) {
    z_left[inode]=par->iz0_all[inode]*dx;
    z_right[inode]=(par->iz0_all[inode]+par->nz_all[inode])*dx;
    z_true_left[inode]=z_left[inode]+dx;
    z_true_right[inode]=z_right[inode]-dx;
    z_bleft_left[inode]=z_left[inode]-dx;
    z_bleft_right[inode]=z_left[inode];
    z_bright_left[inode]=z_right[inode];
    z_bright_right[inode]=z_right[inode]+dx;
    if(inode==0) {
      z_bleft_left[inode]+=par->l_box;
      z_bleft_right[inode]+=par->l_box;
    }
    else if(inode==NNodes-1) {
      z_bright_left[inode]-=par->l_box;
      z_bright_right[inode]-=par->l_box;
    }
  }

  //Figure out which particles need to be moved and which stay
  unsigned long long n_inrange=0,n_inbuffer=0;
  for(ip=0;ip<np_real;ip++) {
    int include=0,send=0;
    if(((z[ip]>=z_left[NodeThis]) && (z[ip]<z_right[NodeThis])) ||
       ((z[ip]>=z_bleft_left[NodeThis]) && (z[ip]<z_bleft_right[NodeThis])) ||
       ((z[ip]>=z_bright_left[NodeThis]) && (z[ip]<z_bright_right[NodeThis])))
      include=1;
    if((z[ip]<=z_true_left[NodeThis]) || (z[ip]>=z_true_right[NodeThis]))
      send=1;

    if(send) {
      if(n_inbuffer>=nbuffer)
	report_error(1,"Out of memory, enlarge buffer\n");
      x_b[n_inbuffer]=x[ip];
      y_b[n_inbuffer]=y[ip];
      z_b[n_inbuffer]=z[ip];
      n_inbuffer++;
    }

    if(include) {
      if(n_inrange>=np_allocated)
	report_error(1,"Out of memory, enlarge buffer\n");
      x[n_inrange]=x[ip];
      y[n_inrange]=y[ip];
      z[n_inrange]=z[ip];
      n_inrange++;
    }
  }

  for(inode=1;inode<NNodes;inode++) {
    int node_to=(NodeThis+inode+NNodes)%NNodes;
    int node_from=(NodeThis-inode+NNodes)%NNodes;
    int i=0,j=n_inbuffer;
    while(1) {
      while(i<n_inbuffer && 
	    !(((z_b[i]>=z_left[node_to]) && (z_b[i]<z_right[node_to])) ||
	      ((z_b[i]>=z_bleft_left[node_to]) && (z_b[i]<z_bleft_right[node_to])) ||
	      ((z_b[i]>=z_bright_left[node_to]) && (z_b[i]<z_bright_right[node_to]))))
	i++;
      while(j>=0 && 
	    (((z_b[j]>=z_left[node_to]) && (z_b[j]<z_right[node_to])) ||
	     ((z_b[j]>=z_bleft_left[node_to]) && (z_b[j]<z_bleft_right[node_to])) ||
	     ((z_b[j]>=z_bright_left[node_to]) && (z_b[j]<z_bright_right[node_to]))))
	j--;

      if(i<j) {
	flouble tmp;
	tmp=x_b[i]; x_b[i]=x_b[j]; x_b[j]=tmp;
	tmp=y_b[i]; y_b[i]=y_b[j]; y_b[j]=tmp;
	tmp=z_b[i]; z_b[i]=z_b[j]; z_b[j]=tmp;
	i++; j--;
      }
      else
	break;
    }

    int nsend,nrecv,nstay;
    if(j==-1) {
      nsend=n_inbuffer;
      nstay=0;
    }
    else if(i==n_inbuffer) {
      nsend=0;
      nstay=n_inbuffer;
    }
    else {
      if(j+1!=i)
	report_error(1,"Error sharing particles\n");
      nstay=i;
      nsend=n_inbuffer-i;
    }

    MPI_Status status;
    flouble *x_send=x_b+nstay;
    flouble *y_send=y_b+nstay;
    flouble *z_send=z_b+nstay;
    int tag=500+5*inode;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Sendrecv(&nsend,1,MPI_INT,node_to,tag,&nrecv,1,MPI_INT,node_from,tag,MPI_COMM_WORLD,&status);
    tag++;

    if(n_inrange+nrecv>np_allocated)
      report_error(1,"Not enough memory, enlarge buffer\n");

    MPI_Sendrecv(x_send,nsend*sizeof(flouble),FLOUBLE_MPI,node_to,tag,
		 x+n_inrange,nrecv*sizeof(flouble),FLOUBLE_MPI,node_from,tag,MPI_COMM_WORLD,&status);
    MPI_Sendrecv(y_send,nsend*sizeof(flouble),FLOUBLE_MPI,node_to,tag,
		 y+n_inrange,nrecv*sizeof(flouble),FLOUBLE_MPI,node_from,tag,MPI_COMM_WORLD,&status);
    MPI_Sendrecv(z_send,nsend*sizeof(flouble),FLOUBLE_MPI,node_to,tag,
		 z+n_inrange,nrecv*sizeof(flouble),FLOUBLE_MPI,node_from,tag,MPI_COMM_WORLD,&status);
    
    n_inrange+=nrecv;
  }

  *np_here=n_inrange;

  free(z_left);
  free(z_right);
  free(z_bleft_left);
  free(z_bleft_right);
  free(z_bright_left);
  free(z_bright_right);
  free(x_b);
  free(y_b);
  free(z_b);

  return;
}

static void lpt_1(ParamCoLoRe *par)
{
  int axis;

  dftw_complex *(cdisp[3]);
  flouble *(disp[3]);
  ptrdiff_t dsize=par->nz_here*((lint)(par->n_grid*(par->n_grid/2+1)));
  ptrdiff_t dsize_buff=(ptrdiff_t)(dsize*(1+par->lpt_buffer_fraction));

  print_info("1-LPT\n");

  print_info(" - Transforming density field\n");
  for(axis=0;axis<3;axis++) {
    cdisp[axis]=dftw_alloc_complex(dsize_buff);
    disp[axis]=(flouble *)cdisp[axis];
  }
  fftw_wrap_r2c(par->n_grid,par->grid_dens,par->grid_dens_f);
  
  print_info(" - Computing displacement field\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par,cdisp)
#endif //_HAVE_OMP
  {
    int ii;
    double dk=2*M_PI/par->l_box;
    double kv[3];

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(ii=0;ii<par->nz_here;ii++) {
      int jj,ii_true;
      ii_true=par->iz0_here+ii;
      if(2*ii_true<=par->n_grid)
	kv[2]=ii_true*dk;
      else
	kv[2]=-(par->n_grid-ii_true)*dk;
      for(jj=0;jj<par->n_grid;jj++) {
	int kk;
	if(2*jj<=par->n_grid)
	  kv[1]=jj*dk;
	else
	  kv[1]=-(par->n_grid-jj)*dk;
	for(kk=0;kk<=par->n_grid/2;kk++) {
	  int ax;
	  double k_mod2;
	  lint index=kk+(par->n_grid/2+1)*((lint)(jj+par->n_grid*ii)); //Grid index for +k
	  if(2*kk<=par->n_grid)
	    kv[0]=kk*dk;
	  else
	    kv[0]=-(par->n_grid-kk)*dk; //This should never happen

	  k_mod2=kv[0]*kv[0]+kv[1]*kv[1]+kv[2]*kv[2];

	  for(ax=2;ax>=0;ax++) {
	    if(k_mod2<=0)
	      cdisp[ax][index]=0;
	    else
	      cdisp[ax][index]=I*kv[ax]*par->grid_dens_f[index]/k_mod2;
	  }
	}
      }
    } //end omp for
  } //end omp parallel

  print_info(" - Transform displacement field\n");
  for(axis=0;axis<3;axis++)
    fftw_wrap_c2r(par->n_grid,cdisp[axis],disp[axis]);

  print_info(" - Computing particle positions\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none) shared(par,disp)
#endif //_HAVE_OMP
  {
    lint iz;
    int ngx=2*(par->n_grid/2+1);
    double dx=par->l_box/par->n_grid;
    double xv[3];

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint indexz=iz*((lint)(ngx*par->n_grid));
      xv[2]=(iz+par->iz0_here+0.0)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint indexy=iy*ngx;
	xv[1]=(iy+0.0)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  int ax;
	  double r,dg;
	  lint index=ix+indexy+indexz;
	  lint index_nopad=ix+par->n_grid*((lint)(iy+par->n_grid*iz));
	  xv[0]=(ix+0.0)*dx-par->pos_obs[0];
	  r=sqrt(xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2]);
	  dg=dgrowth_of_r(par,r);
	  for(ax=0;ax<3;ax++) {
	    flouble p=xv[ax]+dg*disp[ax][index];
	    if(p<0) p+=par->l_box;
	    else if(p>=par->l_box) p-=par->l_box;
	    disp[ax][index_nopad]=p;
	  }
	  par->grid_dens[index]=0;
	}
      }
    } //end omp for
  } //end omp parallel

  unsigned long long np_here;
#ifdef _HAVE_MPI
  print_info(" - Sharing particle positions\n");
  share_particles(par,(unsigned long long)(2*dsize_buff),
		  (unsigned long long)(par->nz_here*((lint)(par->n_grid*par->n_grid))),
		  disp[0],disp[1],disp[2],&np_here);
#else //_HAVE_MPI
  np_nere=par->nz_here*((lint)(par->n_grid*par->n_grid));
#endif //_HAVE_MPI

  print_info(" - Interpolating positions into density field\n");
  pos_2_dens(par,np_here,disp[0],disp[1],disp[2],par->grid_dens);

#ifdef _SPREC
  for(axis=0;axis<3;axis++)
    fftwf_free(cdisp[axis]);
#else //_SPREC
  for(axis=0;axis<3;axis++)
    fftw_free(cdisp[axis]);
#endif //_SPREC

  print_info(" - Normalizing density field\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none) shared(par)
#endif //_HAVE_OMP
  {
    lint iz;
    int ngx=2*(par->n_grid/2+1);
    flouble inv_dens=1.;

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint indexz=iz*((lint)(ngx*par->n_grid));
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint indexy=iy*ngx;
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=ix+indexy+indexz;
	  par->grid_dens[index]=par->grid_dens[index]*inv_dens-1.;
	}
      }
    } //end omp for
  } //end omp parallel

}

//Lognormalizes density grid and puts in lightcone
static void lognormalize(ParamCoLoRe *par)
{
  print_info("Lognormalizin'\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none) shared(par)
#endif //_HAVE_OMP
  {
    lint iz;
    int ngx=2*(par->n_grid/2+1);
    double dx=par->l_box/par->n_grid;

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint indexz=iz*((lint)(ngx*par->n_grid));
      double z0=(iz+par->iz0_here+0.5)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint indexy=iy*ngx;
	double y0=(iy+0.5)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=ix+indexy+indexz;
	  double x0=(ix+0.5)*dx-par->pos_obs[0];
	  double r=sqrt(x0*x0+y0*y0+z0*z0);
	  double dg=dgrowth_of_r(par,r);
	  double delta=par->grid_dens[index];
	  par->grid_dens[index]=exp(dg*(delta-0.5*dg*par->sigma2_gauss))-1;
	}
      }
    }//end omp for
  }//end omp parallel
}

void compute_physical_density_field(ParamCoLoRe *par)
{
  print_info("*** Creating physical matter density\n");
  if(NodeThis==0) timer(0);
  if(par->dens_type==DENS_TYPE_LGNR)
    lognormalize(par);
  else if(par->dens_type==DENS_TYPE_1LPT)
    lpt_1(par);
  else
    report_error(1,"Density type %d not supported\n",par->dens_type);
  if(NodeThis==0) timer(2);
  print_info("\n");
}

static void collect_sigmaz_from_grid(ParamCoLoRe *par,int nz,int idz,
				     double *zarr,double *sarr,unsigned long long *narr)
{
#ifdef _HAVE_OMP
#pragma omp parallel default(none)		\
  shared(par,idz,nz,zarr,narr,sarr)
#endif //_HAVE_OMP
  {
    int iz;
    double dx=par->l_box/par->n_grid;
    double *zarr_thr=my_calloc(nz,sizeof(double));
    double *sarr_thr=my_calloc(nz,sizeof(double));
    unsigned long long *narr_thr=my_calloc(nz,sizeof(unsigned long long));

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      lint iz0=iz*((lint)(2*(par->n_grid/2+1)*par->n_grid));
      double z0=(iz+par->iz0_here+0.5)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	lint iy0=iy*2*(par->n_grid/2+1);
	double y0=(iy+0.5)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  lint index=ix+iy0+iz0;
	  double x0=(ix+0.5)*dx-par->pos_obs[0];
	  double r=sqrt(x0*x0+y0*y0+z0*z0);
	  double redshift=z_of_r(par,r);
	  int ind_z=(int)(redshift*idz)+1;
	  if((ind_z>=0) && (ind_z<nz)) {
	    double d=par->grid_dens[index];
	    narr_thr[iz]++;
	    sarr_thr[iz]+=d*d;
	    zarr_thr[iz]+=redshift;
	  }
	}
      }
    } //end omp for
#ifdef _HAVE_OMP
#pragma omp critical
#endif //_HAVE_OMP
    {
      for(iz=0;iz<nz;iz++) {
	narr[iz]+=narr_thr[iz];
	zarr[iz]+=zarr_thr[iz];
	sarr[iz]+=sarr_thr[iz];
      }
    } //end omp critical

    free(narr_thr);
    free(zarr_thr);
    free(sarr_thr);
  } //end omp parallel
}

static void collect_sigmaz_from_pixels(ParamCoLoRe *par,int nz,int idz,
				       double *zarr,double *sarr,unsigned long long *narr)
{
  int ib;

  for(ib=0;ib<par->n_beams_here;ib++) {
    OnionInfo *oi=par->oi_beams[ib];
#ifdef _HAVE_OMP
#pragma omp parallel default(none) shared(par,oi,ib,zarr,narr,sarr,nz,idz)
#endif //_HAVE_OMP
    {
      int ir;
      double *zarr_thr=my_calloc(nz,sizeof(double));
      double *sarr_thr=my_calloc(nz,sizeof(double));
      unsigned long long *narr_thr=my_calloc(nz,sizeof(unsigned long long));

#ifdef _HAVE_OMP
#pragma omp for nowait schedule(dynamic)
#endif //_HAVE_OMP
      for(ir=0;ir<oi->nr;ir++) {
	int ipix;
	double r=0.5*(oi->rf_arr[ir]+oi->r0_arr[ir]);
	double redshift=z_of_r(par,r);
	int iz=(int)(redshift*idz)+1;

	narr_thr[iz]+=oi->num_pix[ir];
	for(ipix=0;ipix<oi->num_pix[ir];ipix++) {
	  double d=par->dens_beams[ib][ir][ipix];
	  sarr_thr[iz]+=d*d;
	  zarr_thr[iz]+=redshift;
	}
      } //end omp for
#ifdef _HAVE_OMP
#pragma omp critical
#endif //_HAVE_OMP
      {
	for(ir=0;ir<nz;ir++) {
	  narr[ir]+=narr_thr[ir];
	  zarr[ir]+=zarr_thr[ir];
	  sarr[ir]+=sarr_thr[ir];
	}
      } //end omp critical

      free(narr_thr);
      free(zarr_thr);
      free(sarr_thr);
    } //end omp parallel
  }
}

static void collect_sigmaz(ParamCoLoRe *par,int nz,double idz,
			   double *zarr,double *sarr,unsigned long long *narr)
{
  if(par->need_onions)
    collect_sigmaz_from_pixels(par,nz,idz,zarr,sarr,narr);
  else 
    collect_sigmaz_from_grid(par,nz,idz,zarr,sarr,narr);
}

//Computes sigma2(z) for physical density field
void compute_density_normalization(ParamCoLoRe *par)
{
  int nz,iz;
  double idz,df;
  double *zarr,*sarr;
  unsigned long long *narr;
  double zmax=z_of_r(par,par->l_box*0.5);

  print_info("*** Computing normalization of density field\n");
  if(NodeThis==0) timer(0);

  nz=(int)(zmax/DZ_SIGMA)+2;
  idz=(nz-2)/zmax;

  zarr=my_calloc(nz,sizeof(double));
  narr=my_calloc(nz,sizeof(unsigned long long));
  sarr=my_calloc(nz,sizeof(double));

  collect_sigmaz(par,nz,idz,zarr,sarr,narr);

#ifdef _HAVE_MPI
  MPI_Allreduce(MPI_IN_PLACE,narr,nz,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,sarr,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,zarr,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif //_HAVE_MPI

  for(iz=0;iz<nz;iz++) {
    if(narr[iz]>0) {
      zarr[iz]/=narr[iz];
      sarr[iz]/=narr[iz];
    }
  }

  //TODO: only for lognormal
  df=dgrowth_of_r(par,0.5*par->l_box);
  zarr[0]=0;
  sarr[0]=exp(par->sigma2_gauss)-1;
  zarr[nz-1]=z_of_r(par,0.5*par->l_box);
  sarr[nz-1]=exp(df*df*par->sigma2_gauss)-1;

  par->z0_sigma2=zarr[0];
  par->zf_sigma2=zarr[nz-1];
  par->sigma2_0=sarr[0];
  par->sigma2_f=sarr[nz-1];
  par->spline_sigma2_z=gsl_spline_alloc(gsl_interp_cspline,nz);
  par->intacc_sigma2_z=gsl_interp_accel_alloc();
  gsl_spline_init(par->spline_sigma2_z,zarr,sarr,nz);

  if(NodeThis==0) timer(2);
  for(iz=0;iz<nz;iz++) {
    print_info("z=%.3lE, <d^2>=%.3lE, <d^2_s>=%.3lE, n=%d\n",
	       zarr[iz],sqrt(sarr[iz]),sqrt(sigma2_of_z(par,zarr[iz])),narr[iz]);
  }
  print_info("\n");

  free(zarr);
  free(narr);
  free(sarr);
}
