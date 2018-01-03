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
  long ngx=2*(par->n_grid/2+1);

  for(ii=0;ii<np;ii++) {
    long index;
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
      index=i0[0]+ngx*(i0[1]+par->n_grid*i0[2]);
      delta[index]+=1.;
    }
  }
}

static void pos_2_cic(ParamCoLoRe *par,unsigned long long np,
		      flouble *x,flouble *y,flouble *z,flouble *delta)
{
  unsigned long long ii;
  flouble i_agrid=par->n_grid/par->l_box;
  long ngx=2*(par->n_grid/2+1);

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
      delta[i0[0]+ngx*(i0[1]+par->n_grid*i0[2])]+=a0[0]*a0[1]*a0[2];
      delta[i1[0]+ngx*(i0[1]+par->n_grid*i0[2])]+=a1[0]*a0[1]*a0[2];
      delta[i0[0]+ngx*(i1[1]+par->n_grid*i0[2])]+=a0[0]*a1[1]*a0[2];
      delta[i1[0]+ngx*(i1[1]+par->n_grid*i0[2])]+=a1[0]*a1[1]*a0[2];
    }
    if((i1[2]>=0) && (i1[2]<par->nz_here)) {
      delta[i0[0]+ngx*(i0[1]+par->n_grid*i1[2])]+=a0[0]*a0[1]*a1[2];
      delta[i1[0]+ngx*(i0[1]+par->n_grid*i1[2])]+=a1[0]*a0[1]*a1[2];
      delta[i0[0]+ngx*(i1[1]+par->n_grid*i1[2])]+=a0[0]*a1[1]*a1[2];
      delta[i1[0]+ngx*(i1[1]+par->n_grid*i1[2])]+=a1[0]*a1[1]*a1[2];
    }
  }
}

static void pos_2_tsc(ParamCoLoRe *par,unsigned long long np,
		      flouble *x,flouble *y,flouble *z,flouble *delta)
{
  unsigned long long ii;
  flouble i_agrid=par->n_grid/par->l_box;
  long ngx=2*(par->n_grid/2+1);

  for(ii=0;ii<np;ii++) {
    int ax,i0[3],ip[3],im[3];
    flouble a0[3],ap[3],am[3];

    i0[0]=(int)(floorf(x[ii]*i_agrid+0.5));
    i0[1]=(int)(floorf(y[ii]*i_agrid+0.5));
    i0[2]=(int)(floorf(z[ii]*i_agrid+0.5));
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
      delta[im[0]+ngx*(im[1]+par->n_grid*im[2])]+=am[0]*am[1]*am[2];
      delta[i0[0]+ngx*(im[1]+par->n_grid*im[2])]+=a0[0]*am[1]*am[2];
      delta[ip[0]+ngx*(im[1]+par->n_grid*im[2])]+=ap[0]*am[1]*am[2];
      delta[im[0]+ngx*(i0[1]+par->n_grid*im[2])]+=am[0]*a0[1]*am[2];
      delta[i0[0]+ngx*(i0[1]+par->n_grid*im[2])]+=a0[0]*a0[1]*am[2];
      delta[ip[0]+ngx*(i0[1]+par->n_grid*im[2])]+=ap[0]*a0[1]*am[2];
      delta[im[0]+ngx*(ip[1]+par->n_grid*im[2])]+=am[0]*ap[1]*am[2];
      delta[i0[0]+ngx*(ip[1]+par->n_grid*im[2])]+=a0[0]*ap[1]*am[2];
      delta[ip[0]+ngx*(ip[1]+par->n_grid*im[2])]+=ap[0]*ap[1]*am[2];
    }
    if((i0[2]>=0) && (i0[2]<par->nz_here)) {
      delta[im[0]+ngx*(im[1]+par->n_grid*i0[2])]+=am[0]*am[1]*a0[2];
      delta[i0[0]+ngx*(im[1]+par->n_grid*i0[2])]+=a0[0]*am[1]*a0[2];
      delta[ip[0]+ngx*(im[1]+par->n_grid*i0[2])]+=ap[0]*am[1]*a0[2];
      delta[im[0]+ngx*(i0[1]+par->n_grid*i0[2])]+=am[0]*a0[1]*a0[2];
      delta[i0[0]+ngx*(i0[1]+par->n_grid*i0[2])]+=a0[0]*a0[1]*a0[2];
      delta[ip[0]+ngx*(i0[1]+par->n_grid*i0[2])]+=ap[0]*a0[1]*a0[2];
      delta[im[0]+ngx*(ip[1]+par->n_grid*i0[2])]+=am[0]*ap[1]*a0[2];
      delta[i0[0]+ngx*(ip[1]+par->n_grid*i0[2])]+=a0[0]*ap[1]*a0[2];
      delta[ip[0]+ngx*(ip[1]+par->n_grid*i0[2])]+=ap[0]*ap[1]*a0[2];
    }
    if((ip[2]>=0) && (ip[2]<par->nz_here)) {
      delta[im[0]+ngx*(im[1]+par->n_grid*ip[2])]+=am[0]*am[1]*ap[2];
      delta[i0[0]+ngx*(im[1]+par->n_grid*ip[2])]+=a0[0]*am[1]*ap[2];
      delta[ip[0]+ngx*(im[1]+par->n_grid*ip[2])]+=ap[0]*am[1]*ap[2];
      delta[im[0]+ngx*(i0[1]+par->n_grid*ip[2])]+=am[0]*a0[1]*ap[2];
      delta[i0[0]+ngx*(i0[1]+par->n_grid*ip[2])]+=a0[0]*a0[1]*ap[2];
      delta[ip[0]+ngx*(i0[1]+par->n_grid*ip[2])]+=ap[0]*a0[1]*ap[2];
      delta[im[0]+ngx*(ip[1]+par->n_grid*ip[2])]+=am[0]*ap[1]*ap[2];
      delta[i0[0]+ngx*(ip[1]+par->n_grid*ip[2])]+=a0[0]*ap[1]*ap[2];
      delta[ip[0]+ngx*(ip[1]+par->n_grid*ip[2])]+=ap[0]*ap[1]*ap[2];
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

#ifdef _HAVE_MPI
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
						  ((long)(par->n_grid*par->n_grid)));
  flouble *x_b=my_malloc(nbuffer*sizeof(flouble));
  flouble *y_b=my_malloc(nbuffer*sizeof(flouble));
  flouble *z_b=my_malloc(nbuffer*sizeof(flouble));

  for(inode=0;inode<NNodes;inode++) {
    z_left[inode]=par->iz0_all[inode]*dx;
    z_right[inode]=(par->iz0_all[inode]+par->nz_all[inode])*dx;
    z_true_left[inode]=(par->iz0_all[inode]+(par->lpt_interp_type+1)*0.5)*dx;
    z_true_right[inode]=(par->iz0_all[inode]+par->nz_all[inode]-(par->lpt_interp_type+1)*0.5)*dx;
    z_bleft_left[inode]=(par->iz0_all[inode]-(par->lpt_interp_type+1)*0.5)*dx;
    z_bleft_right[inode]=z_left[inode];
    z_bright_left[inode]=z_right[inode];
    z_bright_right[inode]=(par->iz0_all[inode]+par->nz_all[inode]+(par->lpt_interp_type+1)*0.5)*dx;
    if(inode==0) {
      z_bleft_left[inode]=(par->n_grid-(par->lpt_interp_type+1)*0.5)*dx;
      z_bleft_right[inode]=par->n_grid*dx;
    }
    else if(inode==NNodes-1) {
      z_bright_left[inode]=0;
      z_bright_right[inode]=(par->lpt_interp_type+1)*0.5*dx;
    }
#ifdef _DEBUG
    print_info("Node %d: [%lf,%lf] [%lf,%lf] [%lf %lf] [%lf %lf]\n",inode,
	       z_left[inode],z_right[inode],
	       z_true_left[inode],z_true_right[inode],
	       z_bleft_left[inode],z_bleft_right[inode],
	       z_bright_left[inode],z_bright_right[inode]);
#ifdef _HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif //_DEBUG
  }

  //Figure out which particles need to be moved and which stay
  unsigned long long n_inrange=0,n_inbuffer=0;
  for(ip=0;ip<np_real;ip++) {
    int include=0,send=0;
    if(((z[ip]>=z_left[NodeThis]) && (z[ip]<z_right[NodeThis])) ||
       ((z[ip]>=z_bleft_left[NodeThis]) && (z[ip]<z_bleft_right[NodeThis])) ||
       ((z[ip]>=z_bright_left[NodeThis]) && (z[ip]<z_bright_right[NodeThis])))
      include=1;
    if((z[ip]<z_true_left[NodeThis]) || (z[ip]>=z_true_right[NodeThis]))
      send=1;

    if(send) {
      if(n_inbuffer>=nbuffer)
	report_error(1,"Out of memory, enlarge buffer %llu (%llu,%llu) %llu\n",
		     np_allocated,n_inbuffer,nbuffer,n_inrange);
      x_b[n_inbuffer]=x[ip];
      y_b[n_inbuffer]=y[ip];
      z_b[n_inbuffer]=z[ip];
      n_inbuffer++;
    }

    if(include) {
      if(n_inrange>=np_allocated)
	report_error(1,"Out of memory, enlarge buffer %llu %llu %llu\n",np_allocated,n_inbuffer,n_inrange);
      x[n_inrange]=x[ip];
      y[n_inrange]=y[ip];
      z[n_inrange]=z[ip];
      n_inrange++;
    }
  }

#ifdef _DEBUG
  printf("Node %d: %llu %llu %llu %llu\n",NodeThis,np_allocated,n_inrange,nbuffer,n_inbuffer);
#ifdef _HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif //_DEBUG

  for(inode=1;inode<NNodes;inode++) {
    int node_to=(NodeThis+inode+NNodes)%NNodes;
    int node_from=(NodeThis-inode+NNodes)%NNodes;
    int i=0,j=n_inbuffer-1;

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

#ifdef _HAVE_MPI
    flouble *x_send=x_b+nstay;
    flouble *y_send=y_b+nstay;
    flouble *z_send=z_b+nstay;
#endif //_HAVE_MPI
    int tag=500+5*inode;

#ifdef _HAVE_MPI
    MPI_Status status;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Sendrecv(&nsend,1,MPI_INT,node_to,tag,&nrecv,1,MPI_INT,node_from,tag,MPI_COMM_WORLD,&status);
#endif //_HAVE_MPI
    tag++;
    

#ifdef _DEBUG
    printf("%d. Node %d: send %d particles to node %d, receive %d particles from node %d\n",
	   inode,NodeThis,nsend,node_to,nrecv,node_from);
#ifdef _HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif //_HAVE_MPI
#endif //_DEBUG

    if(n_inrange+nrecv>np_allocated)
      report_error(1,"Not enough memory, enlarge buffer\n");

#ifdef _HAVE_MPI
    MPI_Sendrecv(x_send,nsend*sizeof(flouble),FLOUBLE_MPI,node_to,tag,
		 x+n_inrange,nrecv*sizeof(flouble),FLOUBLE_MPI,node_from,tag,MPI_COMM_WORLD,&status);
    MPI_Sendrecv(y_send,nsend*sizeof(flouble),FLOUBLE_MPI,node_to,tag,
		 y+n_inrange,nrecv*sizeof(flouble),FLOUBLE_MPI,node_from,tag,MPI_COMM_WORLD,&status);
    MPI_Sendrecv(z_send,nsend*sizeof(flouble),FLOUBLE_MPI,node_to,tag,
		 z+n_inrange,nrecv*sizeof(flouble),FLOUBLE_MPI,node_from,tag,MPI_COMM_WORLD,&status);
#endif //_HAVE_MPI

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
#endif //_HAVE_MPI

static void lpt_1(ParamCoLoRe *par)
{
  int axis;

  dftw_complex *(cdisp[3]);
  flouble *(disp[3]);
  ptrdiff_t dsize=par->nz_here*((long)(par->n_grid*(par->n_grid/2+1)));
  ptrdiff_t dsize_buff=(ptrdiff_t)(dsize*(1+par->lpt_buffer_fraction));

  print_info(" 1LPT\n");

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
    double fftnorm=(double)(par->n_grid*((long)(par->n_grid*par->n_grid)));

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
	  long index=kk+(par->n_grid/2+1)*((long)(jj+par->n_grid*ii)); //Grid index for +k
	  if(2*kk<=par->n_grid)
	    kv[0]=kk*dk;
	  else
	    kv[0]=-(par->n_grid-kk)*dk; //This should never happen

	  k_mod2=fftnorm*(kv[0]*kv[0]+kv[1]*kv[1]+kv[2]*kv[2]);

	  for(ax=0;ax<3;ax++) {
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
#ifdef _DEBUG
  double d_sigma2_1=0;
  double d_mean_1[3]={0,0,0};
#endif //_DEBUG
#ifdef _HAVE_OMP
#ifdef _DEBUG
#pragma omp parallel default(none) shared(par,disp,d_mean_1,d_sigma2_1)
#else //_DEBUG
#pragma omp parallel default(none) shared(par,disp)
#endif //_DEBUG
#endif //_HAVE_OMP
  {
    long iz;
    int ngx=2*(par->n_grid/2+1);
    flouble dx=par->l_box/par->n_grid;
    flouble xv[3];

#ifdef _DEBUG
    double d_sigma2_1_thr=0;
    double d_mean_1_thr[3]={0,0,0};
#endif //_DEBUG

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      xv[2]=(iz+par->iz0_here+0.0)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	xv[1]=(iy+0.0)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  int ax;
	  double r,dg;
	  long index=ix+indexy+indexz;
	  xv[0]=(ix+0.0)*dx-par->pos_obs[0];
	  r=sqrt(xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2]);
	  dg=get_bg(par,r,BG_D1,0);
	  for(ax=0;ax<3;ax++) {
#ifdef _DEBUG
	    d_mean_1_thr[ax]+=disp[ax][index];
	    d_sigma2_1_thr+=disp[ax][index]*disp[ax][index];
#endif //_DEBUG
	    flouble p=xv[ax]+dg*disp[ax][index]+par->pos_obs[ax];
	    if(p<0) p+=par->l_box;
	    if(p>=par->l_box) p-=par->l_box;
	    disp[ax][index]=p;
	  }
	  par->grid_dens[index]=0;
	}
      }
    } //end omp for
#ifdef _DEBUG
#ifdef _HAVE_OMP
#pragma omp critical
#endif //_HAVE_OMP
    {
      d_mean_1[0]+=d_mean_1_thr[0];
      d_mean_1[1]+=d_mean_1_thr[1];
      d_mean_1[2]+=d_mean_1_thr[2];
      d_sigma2_1+=d_sigma2_1_thr;
    }
#endif //_DEBUG
  } //end omp parallel

#ifdef _DEBUG
#ifdef _HAVE_MPI
  if(NodeThis==0) 
    MPI_Reduce(MPI_IN_PLACE,d_mean_1,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  else 
    MPI_Reduce(d_mean_1,d_mean_1,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(NodeThis==0) 
    MPI_Reduce(MPI_IN_PLACE,&d_sigma2_1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  else 
    MPI_Reduce(&d_sigma2_1,&d_sigma2_1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#endif //_HAVE_MPI
  if(NodeThis==0) {
    d_mean_1[0]/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
    d_mean_1[1]/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
    d_mean_1[2]/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
    d_sigma2_1/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
  
    print_info(" 1st-order displacement : [%lE,%lE,%lE] %lE\n",
	       d_mean_1[0],d_mean_1[1],d_mean_1[2],sqrt(d_sigma2_1));
  }
#endif //_DEBUG

  print_info(" - Undoing padding\n");
  {
    long iz;
    int ngx=2*(par->n_grid/2+1);

    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      long indexz_unpadded=iz*((long)(par->n_grid*par->n_grid));
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	long indexy_unpadded=iy*par->n_grid;
	for(ix=0;ix<par->n_grid;ix++) {
	  int ax;
	  long index=ix+indexy+indexz;
	  long index_unpadded=ix+indexy_unpadded+indexz_unpadded;
	  for(ax=0;ax<3;ax++)
	    disp[ax][index_unpadded]=disp[ax][index];
	}
      }
    }
  }

  unsigned long long np_here;
  if(par->output_lpt) {
    np_here=par->nz_here*((long)(par->n_grid*par->n_grid));
    print_info(" - Writing LPT positions\n");
    write_lpt(par,np_here,disp[0],disp[1],disp[2]);
  }

#ifdef _HAVE_MPI
  print_info(" - Sharing particle positions\n");
  share_particles(par,(unsigned long long)(2*dsize_buff),
		  (unsigned long long)(par->nz_here*((long)(par->n_grid*par->n_grid))),
		  disp[0],disp[1],disp[2],&np_here);
#else //_HAVE_MPI
  np_here=par->nz_here*((long)(par->n_grid*par->n_grid));
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
#ifdef _DEBUG
  double numtot=0;
#endif //_DEBUG
#ifdef _HAVE_OMP
#ifdef _DEBUG
#pragma omp parallel default(none) shared(par,numtot)
#else //_DEBUG
#pragma omp parallel default(none) shared(par)
#endif //_DEBUG
#endif //_HAVE_OMP
  {
    long iz;
    int ngx=2*(par->n_grid/2+1);
    flouble inv_dens=1.;
#ifdef _DEBUG
    double numtot_thr=0;
#endif //_DEBUG

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	for(ix=0;ix<par->n_grid;ix++) {
	  long index=ix+indexy+indexz;
#ifdef _DEBUG
	  numtot_thr+=par->grid_dens[index];
#endif //_DEBUG
	  par->grid_dens[index]=par->grid_dens[index]*inv_dens-1.;
	}
      }
    } //end omp for
#ifdef _DEBUG
#ifdef _HAVE_OMP
#pragma omp critical
#endif //_HAVE_OMP
    {
      numtot+=numtot_thr;
    }
#endif //_DEBUG
  } //end omp parallel

#ifdef _DEBUG
#ifdef _HAVE_MPI
  if(NodeThis==0)
    MPI_Reduce(MPI_IN_PLACE,&numtot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  else
    MPI_Reduce(&numtot,&numtot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#endif //_HAVE_MPI
  print_info(" Total density : %lE\n",numtot-(double)(par->n_grid*((long)(par->n_grid*par->n_grid))));
#endif //_DEBUG
}

static void lpt_2(ParamCoLoRe *par)
{
  int axis;

  dftw_complex *(cdisp[3]),*(cdigrad[6]);
  flouble *(disp[3]),*(digrad[6]);
  ptrdiff_t dsize=par->nz_here*((long)(par->n_grid*(par->n_grid/2+1)));
  ptrdiff_t dsize_buff=(ptrdiff_t)(dsize*(1+par->lpt_buffer_fraction));

  print_info(" 2LPT\n");

  print_info(" - Transforming density field\n");
  for(axis=0;axis<2;axis++) {
    cdisp[axis]=dftw_alloc_complex(dsize_buff);
    disp[axis]=(flouble *)cdisp[axis];
  }
  for(axis=0;axis<6;axis++) {
    cdigrad[axis]=dftw_alloc_complex(dsize_buff);
    digrad[axis]=(flouble *)cdigrad[axis];
  }
  cdisp[2]=par->grid_dens_f;
  disp[2]=par->grid_dens;
  fftw_wrap_r2c(par->n_grid,par->grid_dens,par->grid_dens_f);

  print_info(" - Computing displacement field\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par,cdisp,cdigrad)
#endif //_HAVE_OMP
  {
    int ii;
    double dk=2*M_PI/par->l_box;
    double kv[3];
    double fftnorm=(double)(par->n_grid*((long)(par->n_grid*par->n_grid)));

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
	  long index=kk+(par->n_grid/2+1)*((long)(jj+par->n_grid*ii)); //Grid index for +k
	  if(2*kk<=par->n_grid)
	    kv[0]=kk*dk;
	  else
	    kv[0]=-(par->n_grid-kk)*dk; //This should never happen

	  k_mod2=fftnorm*(kv[0]*kv[0]+kv[1]*kv[1]+kv[2]*kv[2]);

	  for(ax=0;ax<3;ax++) {
	    if(k_mod2<=0)
	      cdisp[ax][index]=0;
	    else
	      cdisp[ax][index]=I*kv[ax]*par->grid_dens_f[index]/k_mod2;
	  }

	  cdigrad[0][index]=I*kv[0]*cdisp[0][index]; //SIGN
	  cdigrad[1][index]=I*kv[1]*cdisp[0][index];
	  cdigrad[2][index]=I*kv[2]*cdisp[0][index];
	  cdigrad[3][index]=I*kv[1]*cdisp[1][index];
	  cdigrad[4][index]=I*kv[2]*cdisp[1][index];
	  cdigrad[5][index]=I*kv[2]*cdisp[2][index];
	}
      }
    } //end omp for
  } //end omp parallel

  print_info(" - Transform digradient\n");
  for(axis=0;axis<6;axis++)
    fftw_wrap_c2r(par->n_grid,cdigrad[axis],digrad[axis]);

  print_info(" - Computing second-order potential\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par,digrad)
#endif //_HAVE_OMP
  {
    int iz;
    int ngx=2*(par->n_grid/2+1);

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	for(ix=0;ix<par->n_grid;ix++) {
	  long index=ix+indexy+indexz; //SIGN
	  digrad[5][index]=
	    digrad[0][index]*digrad[3][index]+//xx*yy
	    digrad[0][index]*digrad[5][index]+//xx*zz
	    digrad[3][index]*digrad[5][index]-//+yy*zz
	    digrad[1][index]*digrad[1][index]-//-xy*xy
	    digrad[2][index]*digrad[2][index]-//-xz*xz
	    digrad[4][index]*digrad[4][index];//-yz*yz
	}
      }
    } //end omp for
  } //end omp parallel

  print_info(" - Transforming second-order potential\n");
  fftw_wrap_r2c(par->n_grid,digrad[5],cdigrad[5]);

  print_info(" - Computing 2nd-order displacement field\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none) \
  shared(par,cdigrad)
#endif //_HAVE_OMP
  {
    int ii;
    double dk=2*M_PI/par->l_box;
    double kv[3];
    double fftnorm=(double)(par->n_grid*((long)(par->n_grid*par->n_grid)));

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
	  long index=kk+(par->n_grid/2+1)*((long)(jj+par->n_grid*ii)); //Grid index for +k
	  if(2*kk<=par->n_grid)
	    kv[0]=kk*dk;
	  else
	    kv[0]=-(par->n_grid-kk)*dk; //This should never happen

	  k_mod2=fftnorm*(kv[0]*kv[0]+kv[1]*kv[1]+kv[2]*kv[2]);

	  for(ax=0;ax<3;ax++) {
	    if(k_mod2<=0)
	      cdigrad[ax][index]=0;
	    else
	      cdigrad[ax][index]=-I*kv[ax]*cdigrad[5][index]/k_mod2; //SIGN, 1/k^2, normalization
	  }
	}
      }
    } //end omp for
  } //end omp parallel

  print_info(" - Transform 1st- and 2nd-order displacement fields\n");
  for(axis=0;axis<3;axis++) {
    fftw_wrap_c2r(par->n_grid,cdisp[axis],disp[axis]);
    fftw_wrap_c2r(par->n_grid,cdigrad[axis],digrad[axis]);
  }

  print_info(" - Computing particle positions\n");
#ifdef _DEBUG
  double d_sigma2_1=0;
  double d_mean_1[3]={0,0,0};
  double d_sigma2_2=0;
  double d_mean_2[3]={0,0,0};
#endif //_DEBUG
#ifdef _HAVE_OMP
#ifdef _DEBUG
#pragma omp parallel default(none) shared(par,disp,digrad,d_mean_1,d_mean_2,d_sigma2_1,d_sigma2_2)
#else //_DEBUG
#pragma omp parallel default(none) shared(par,disp,digrad)
#endif //_DEBUG
#endif //_HAVE_OMP
  {
    long iz;
    int ngx=2*(par->n_grid/2+1);
    flouble dx=par->l_box/par->n_grid;
    flouble xv[3];

#ifdef _DEBUG
    double d_sigma2_1_thr=0;
    double d_mean_1_thr[3]={0,0,0};
    double d_sigma2_2_thr=0;
    double d_mean_2_thr[3]={0,0,0};
#endif //_DEBUG

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      xv[2]=(iz+par->iz0_here+0.0)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	xv[1]=(iy+0.0)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  int ax;
	  double r,dg,d2g;
	  long index=ix+indexy+indexz;
	  long index_nopad=ix+par->n_grid*((long)(iy+par->n_grid*iz));
	  xv[0]=(ix+0.0)*dx-par->pos_obs[0];
	  r=sqrt(xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2]);
	  dg=get_bg(par,r,BG_D1,0);
	  d2g=get_bg(par,r,BG_D2,0);
	  for(ax=0;ax<3;ax++) {
#ifdef _DEBUG
	    d_mean_1_thr[ax]+=disp[ax][index];
	    d_mean_2_thr[ax]+=digrad[ax][index];
	    d_sigma2_1_thr+=disp[ax][index]*disp[ax][index];
	    d_sigma2_2_thr+=digrad[ax][index]*digrad[ax][index];
#endif //_DEBUG
	    flouble p=xv[ax]+dg*disp[ax][index]+d2g*digrad[ax][index]+par->pos_obs[ax];
	    if(p<0) p+=par->l_box;
	    if(p>=par->l_box) p-=par->l_box;
	    digrad[3+ax][index_nopad]=p;
	  }
	  par->grid_dens[index]=0;
	}
      }
    } //end omp for
#ifdef _DEBUG
#ifdef _HAVE_OMP
#pragma omp critical
#endif //_HAVE_OMP
    {
      d_mean_1[0]+=d_mean_1_thr[0];
      d_mean_1[1]+=d_mean_1_thr[1];
      d_mean_1[2]+=d_mean_1_thr[2];
      d_mean_2[0]+=d_mean_2_thr[0];
      d_mean_2[1]+=d_mean_2_thr[1];
      d_mean_2[2]+=d_mean_2_thr[2];
      d_sigma2_1+=d_sigma2_1_thr;
      d_sigma2_2+=d_sigma2_2_thr;
    }
#endif //_DEBUG
  } //end omp parallel

#ifdef _DEBUG
#ifdef _HAVE_MPI
  if(NodeThis==0) 
    MPI_Reduce(MPI_IN_PLACE,d_mean_1,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  else 
    MPI_Reduce(d_mean_1,d_mean_1,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(NodeThis==0) 
    MPI_Reduce(MPI_IN_PLACE,d_mean_2,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  else 
    MPI_Reduce(d_mean_2,d_mean_2,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(NodeThis==0) 
    MPI_Reduce(MPI_IN_PLACE,&d_sigma2_1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  else 
    MPI_Reduce(&d_sigma2_1,&d_sigma2_1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(NodeThis==0) 
    MPI_Reduce(MPI_IN_PLACE,&d_sigma2_2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  else 
    MPI_Reduce(&d_sigma2_2,&d_sigma2_2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#endif //_HAVE_MPI
  if(NodeThis==0) {
    d_mean_1[0]/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
    d_mean_1[1]/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
    d_mean_1[2]/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
    d_mean_2[0]/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
    d_mean_2[1]/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
    d_mean_2[2]/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
    d_sigma2_1/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
    d_sigma2_2/=(par->n_grid*((long)(par->n_grid*par->n_grid)));
  
    print_info(" 1st-order displacement : [%lE,%lE,%lE] %lE\n",
	       d_mean_1[0],d_mean_1[1],d_mean_1[2],sqrt(d_sigma2_1));
    print_info(" 2nd-order displacement : [%lE,%lE,%lE] %lE\n",
	       d_mean_2[0],d_mean_2[1],d_mean_2[2],sqrt(d_sigma2_2));
  }
#endif //_DEBUG

#ifdef _SPREC
  for(axis=0;axis<2;axis++)
    fftwf_free(cdisp[axis]);
  for(axis=0;axis<3;axis++)
    fftwf_free(cdigrad[axis]);
#else //_SPREC
  for(axis=0;axis<2;axis++)
    fftw_free(cdisp[axis]);
  for(axis=0;axis<3;axis++)
    fftw_free(cdigrad[axis]);
#endif //_SPREC

  unsigned long long np_here;
  if(par->output_lpt) {
    np_here=par->nz_here*((long)(par->n_grid*par->n_grid));
    print_info(" - Writing LPT positions\n");
    write_lpt(par,np_here,digrad[3],digrad[4],digrad[5]);
  }

#ifdef _HAVE_MPI
  print_info(" - Sharing particle positions\n");
  share_particles(par,(unsigned long long)(2*dsize_buff),
		  (unsigned long long)(par->nz_here*((long)(par->n_grid*par->n_grid))),
		  digrad[3],digrad[4],digrad[5],&np_here);
#else //_HAVE_MPI
  np_here=par->nz_here*((long)(par->n_grid*par->n_grid));
#endif //_HAVE_MPI

  print_info(" - Interpolating positions into density field\n");
  pos_2_dens(par,np_here,digrad[3],digrad[4],digrad[5],par->grid_dens);

#ifdef _SPREC
  for(axis=0;axis<3;axis++)
    fftwf_free(cdigrad[3+axis]);
#else //_SPREC
  for(axis=0;axis<3;axis++)
    fftw_free(cdigrad[3+axis]);
#endif //_SPREC

  print_info(" - Normalizing density field\n");
#ifdef _DEBUG
  double numtot=0;
#endif //_DEBUG
#ifdef _HAVE_OMP
#ifdef _DEBUG
#pragma omp parallel default(none) shared(par,numtot)
#else //_DEBUG
#pragma omp parallel default(none) shared(par)
#endif //_DEBUG
#endif //_HAVE_OMP
  {
    long iz;
    int ngx=2*(par->n_grid/2+1);
    flouble inv_dens=1.;
#ifdef _DEBUG
    double numtot_thr=0;
#endif //_DEBUG

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	for(ix=0;ix<par->n_grid;ix++) {
	  long index=ix+indexy+indexz;
#ifdef _DEBUG
	  numtot_thr+=par->grid_dens[index];
#endif //_DEBUG
	  par->grid_dens[index]=par->grid_dens[index]*inv_dens-1.;
	}
      }
    } //end omp for
#ifdef _DEBUG
#ifdef _HAVE_OMP
#pragma omp critical
#endif //_HAVE_OMP
    {
      numtot+=numtot_thr;
    }
#endif //_DEBUG
  } //end omp parallel

#ifdef _DEBUG
#ifdef _HAVE_MPI
  if(NodeThis==0)
    MPI_Reduce(MPI_IN_PLACE,&numtot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  else
    MPI_Reduce(&numtot,&numtot,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#endif //_HAVE_MPI
  print_info(" Total density : %lE\n",numtot-(double)(par->n_grid*((long)(par->n_grid*par->n_grid))));
#endif //_DEBUG
}

//Lognormalizes density grid and puts in lightcone
static void densclip(ParamCoLoRe *par)
{
  print_info("Lognormalizin'\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none) shared(par)
#endif //_HAVE_OMP
  {
    long iz;
    int ngx=2*(par->n_grid/2+1);
    flouble dx=par->l_box/par->n_grid;

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      flouble z0=(iz+par->iz0_here+0.5)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	flouble y0=(iy+0.5)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  long index=ix+indexy+indexz;
	  flouble x0=(ix+0.5)*dx-par->pos_obs[0];
	  double r=sqrt(x0*x0+y0*y0+z0*z0);
	  double dg=get_bg(par,r,BG_D1,0);
	  double delta=par->grid_dens[index];
	  par->grid_dens[index]=fmax(1+dg*delta,0)-1;
	}
      }
    }//end omp for
  }//end omp parallel
}

//Lognormalizes density grid and puts in lightcone
static void lognormalize(ParamCoLoRe *par)
{
  print_info("Lognormalizin'\n");
#ifdef _HAVE_OMP
#pragma omp parallel default(none) shared(par)
#endif //_HAVE_OMP
  {
    long iz;
    int ngx=2*(par->n_grid/2+1);
    flouble dx=par->l_box/par->n_grid;

#ifdef _HAVE_OMP
#pragma omp for schedule(static)
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long indexz=iz*((long)(ngx*par->n_grid));
      flouble z0=(iz+par->iz0_here+0.5)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long indexy=iy*ngx;
	flouble y0=(iy+0.5)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  long index=ix+indexy+indexz;
	  flouble x0=(ix+0.5)*dx-par->pos_obs[0];
	  double r=sqrt(x0*x0+y0*y0+z0*z0);
	  double dg=get_bg(par,r,BG_D1,0);
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
  else if(par->dens_type==DENS_TYPE_2LPT)
    lpt_2(par);
  else if(par->dens_type==DENS_TYPE_CLIP)
    densclip(par);
  else
    report_error(1,"Density type %d not supported\n",par->dens_type);

  if(NodeThis==0) timer(2);
  print_info("\n");

  if(par->output_density)
    write_density_grid(par,"lightcone");
}

static void collect_density_normalization_from_grid(ParamCoLoRe *par,int nz,double idz,double *zarr,
						    unsigned long long *narr,
						    double **norm_srcs_arr,double **norm_imap_arr)
{
#ifdef _HAVE_OMP
#pragma omp parallel default(none)				\
  shared(par,idz,nz,zarr,narr,norm_srcs_arr,norm_imap_arr)
#endif //_HAVE_OMP
  {
    int iz,ipop;
    flouble dx=par->l_box/par->n_grid;
    double *zarr_thr=my_calloc(nz,sizeof(double));
    unsigned long long *narr_thr=my_calloc(nz,sizeof(unsigned long long));
    double **norm_srcs_arr_thr=my_malloc(par->n_srcs*sizeof(double *));
    double **norm_imap_arr_thr=my_malloc(par->n_imap*sizeof(double *));
    for(ipop=0;ipop<par->n_srcs;ipop++)
      norm_srcs_arr_thr[ipop]=my_calloc(nz,sizeof(double));
    for(ipop=0;ipop<par->n_imap;ipop++)
      norm_imap_arr_thr[ipop]=my_calloc(nz,sizeof(double));

#ifdef _HAVE_OMP
#pragma omp for
#endif //_HAVE_OMP
    for(iz=0;iz<par->nz_here;iz++) {
      int iy;
      long iz0=iz*((long)(2*(par->n_grid/2+1)*par->n_grid));
      flouble z0=(iz+par->iz0_here+0.5)*dx-par->pos_obs[2];
      for(iy=0;iy<par->n_grid;iy++) {
	int ix;
	long iy0=iy*2*(par->n_grid/2+1);
	flouble y0=(iy+0.5)*dx-par->pos_obs[1];
	for(ix=0;ix<par->n_grid;ix++) {
	  long index=ix+iy0+iz0;
	  flouble x0=(ix+0.5)*dx-par->pos_obs[0];
	  double r=sqrt(x0*x0+y0*y0+z0*z0);
	  double redshift=get_bg(par,r,BG_Z,0);
	  int ind_z=(int)(redshift*idz)+1;
	  if((ind_z>=0) && (ind_z<nz)) {
	    double d=par->grid_dens[index];
	    narr_thr[ind_z]++;
	    zarr_thr[ind_z]+=redshift;
	    for(ipop=0;ipop<par->n_srcs;ipop++)
	      norm_srcs_arr_thr[ipop][ind_z]+=bias_model(d,get_bg(par,r,BG_BZ_SRCS,ipop));
	    for(ipop=0;ipop<par->n_imap;ipop++)
	      norm_imap_arr_thr[ipop][ind_z]+=bias_model(d,get_bg(par,r,BG_BZ_IMAP,ipop));
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
	for(ipop=0;ipop<par->n_srcs;ipop++)
	  norm_srcs_arr[ipop][iz]+=norm_srcs_arr_thr[ipop][iz];
	for(ipop=0;ipop<par->n_imap;ipop++)
	  norm_imap_arr[ipop][iz]+=norm_imap_arr_thr[ipop][iz];
      }
    } //end omp critical

    free(narr_thr);
    free(zarr_thr);
    for(ipop=0;ipop<par->n_srcs;ipop++)
      free(norm_srcs_arr_thr[ipop]);
    free(norm_srcs_arr_thr);
    for(ipop=0;ipop<par->n_imap;ipop++)
      free(norm_imap_arr_thr[ipop]);
    free(norm_imap_arr_thr);
  } //end omp parallel
}

static void collect_density_normalization(ParamCoLoRe *par,int nz,double idz,double *zarr,
					  unsigned long long *narr,
					  double **norm_srcs_arr,double **norm_imap_arr)
{
  //  if(par->need_onions)
  //    collect_density_normalization_from_pixels(par,nz,idz,zarr,narr,norm_srcs_arr,norm_imap_arr);
  //  else 
  //    collect_density_normalization_from_grid(par,nz,idz,zarr,narr,norm_srcs_arr,norm_imap_arr);
  collect_density_normalization_from_grid(par,nz,idz,zarr,narr,norm_srcs_arr,norm_imap_arr);
}

//Computes sigma2(z) for physical density field
void compute_density_normalization(ParamCoLoRe *par)
{
  int nz,iz,ii,ipop;
  double idz;
  double *zarr,**norm_imap_arr,**norm_srcs_arr;
  unsigned long long *narr;
  double zmax=get_bg(par,par->l_box*0.5,BG_Z,0);
  gsl_spline *spline_norm_srcs[NPOP_MAX];
  gsl_interp_accel *intacc_srcs=gsl_interp_accel_alloc(); 
  gsl_spline *spline_norm_imap[NPOP_MAX];
  gsl_interp_accel *intacc_imap=gsl_interp_accel_alloc();

  print_info("*** Computing normalization of density field\n");
  if(NodeThis==0) timer(0);

  nz=(int)(zmax/DZ_SIGMA)+2;
  idz=(nz-2)/zmax;

  zarr=my_calloc(nz,sizeof(double));
  narr=my_calloc(nz,sizeof(unsigned long long));
  norm_srcs_arr=my_malloc(par->n_srcs*sizeof(double *));
  norm_imap_arr=my_malloc(par->n_imap*sizeof(double *));
  for(ipop=0;ipop<par->n_srcs;ipop++)
    norm_srcs_arr[ipop]=my_calloc(nz,sizeof(double));
  for(ipop=0;ipop<par->n_imap;ipop++)
    norm_imap_arr[ipop]=my_calloc(nz,sizeof(double));

  collect_density_normalization(par,nz,idz,zarr,narr,norm_srcs_arr,norm_imap_arr);

#ifdef _HAVE_MPI
  MPI_Allreduce(MPI_IN_PLACE,narr,nz,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  for(ipop=0;ipop<par->n_srcs;ipop++)
    MPI_Allreduce(MPI_IN_PLACE,norm_srcs_arr[ipop],nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(ipop=0;ipop<par->n_imap;ipop++)
    MPI_Allreduce(MPI_IN_PLACE,norm_imap_arr[ipop],nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,zarr,nz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif //_HAVE_MPI

  for(iz=0;iz<nz;iz++) {
    if(narr[iz]>0) {
      zarr[iz]/=narr[iz];
      for(ipop=0;ipop<par->n_srcs;ipop++)
	norm_srcs_arr[ipop][iz]=narr[iz]/norm_srcs_arr[ipop][iz];
      for(ipop=0;ipop<par->n_imap;ipop++)
	norm_imap_arr[ipop][iz]=narr[iz]/norm_imap_arr[ipop][iz];
    }
  }

  zarr[0]=0;
  zarr[nz-1]=get_bg(par,0.5*par->l_box,BG_Z,0);
  for(ipop=0;ipop<par->n_srcs;ipop++) {
    norm_srcs_arr[ipop][0]=norm_srcs_arr[ipop][1];
    norm_srcs_arr[ipop][nz-1]=norm_srcs_arr[ipop][nz-2];
  }
  for(ipop=0;ipop<par->n_imap;ipop++) {
    norm_imap_arr[ipop][0]=norm_imap_arr[ipop][1];
    norm_imap_arr[ipop][nz-1]=norm_imap_arr[ipop][nz-2];
  }

  par->z0_norm=zarr[0];
  par->zf_norm=zarr[nz-1];
  for(ipop=0;ipop<par->n_srcs;ipop++) {
    par->norm_srcs_0[ipop]=norm_srcs_arr[ipop][0];
    par->norm_srcs_f[ipop]=norm_srcs_arr[ipop][nz-1];
    spline_norm_srcs[ipop]=gsl_spline_alloc(gsl_interp_linear,nz);
    gsl_spline_init(spline_norm_srcs[ipop],zarr,norm_srcs_arr[ipop],nz);
    par->srcs_norm_arr[ipop]=my_malloc(NA*sizeof(double));
  }
  for(ipop=0;ipop<par->n_imap;ipop++) {
    par->norm_imap_0[ipop]=norm_imap_arr[ipop][0];
    par->norm_imap_f[ipop]=norm_imap_arr[ipop][nz-1];
    spline_norm_imap[ipop]=gsl_spline_alloc(gsl_interp_linear,nz);
    gsl_spline_init(spline_norm_imap[ipop],zarr,norm_imap_arr[ipop],nz);
    par->imap_norm_arr[ipop]=my_malloc(NA*sizeof(double));
  }

  for(ii=0;ii<NA;ii++) {
    double z=get_bg(par,par->r_arr_r2z[ii],BG_Z,0);
    for(ipop=0;ipop<par->n_srcs;ipop++) {
      double nm;
      if(z<par->z0_norm)
	nm=par->norm_srcs_0[ipop];
      else if(z>=par->zf_norm)
	nm=par->norm_srcs_f[ipop];
      else
	nm=gsl_spline_eval(spline_norm_srcs[ipop],z,intacc_srcs);
      par->srcs_norm_arr[ipop][ii]=nm;
    }
    for(ipop=0;ipop<par->n_imap;ipop++) {
      double nm;
      if(z<par->z0_norm)
	nm=par->norm_imap_0[ipop];
      else if(z>=par->zf_norm)
	nm=par->norm_imap_f[ipop];
      else
	nm=gsl_spline_eval(spline_norm_imap[ipop],z,intacc_imap);
      par->imap_norm_arr[ipop][ii]=nm;
    }
  }

#ifdef _DEBUG
  for(iz=0;iz<nz;iz++) {
    double rz=r_of_z(par,zarr[iz]);
    print_info("z=%.3lE, ",zarr[iz]);
    for(ipop=0;ipop<par->n_srcs;ipop++)
      print_info("<d^2_%d>=%.3lE, ",ipop,get_bg(par,rz,BG_NORM_SRCS,ipop));
    for(ipop=0;ipop<par->n_imap;ipop++)
      print_info("<d^2_%d>=%.3lE, ",ipop,get_bg(par,rz,BG_NORM_IMAP,ipop));
    print_info("%011llu\n",narr[iz]);
  }
#endif //_DEBUG
  if(NodeThis==0) timer(2);
  print_info("\n");

  free(zarr);
  free(narr);
  for(ipop=0;ipop<par->n_srcs;ipop++) {
    free(norm_srcs_arr[ipop]);
    gsl_spline_free(spline_norm_srcs[ipop]);
  }
  gsl_interp_accel_free(intacc_srcs);
  for(ipop=0;ipop<par->n_imap;ipop++) {
    free(norm_imap_arr[ipop]);
    gsl_spline_free(spline_norm_imap[ipop]);
  }
  gsl_interp_accel_free(intacc_imap);
}
