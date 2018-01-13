///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso, Anze Slosar                        //
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
#include "fftlog.h"

void write_predictions(ParamCoLoRe *par)
{
  if ((!par->do_srcs) && (!par->do_imap)) return;
  if (NodeThis!=0) return;
  print_info("*** Writing predictions\n");
  // first generate k array, sufficiently finely spaced
  // note that we need to sufficiently pad on both ends
  const int Nk=10000;
  const double kmin=1e-3;
  const double kmax=50;;
  const double kminout=kmin;
  const double kmaxout=kmax;
  const double rminout=0.5;
  const double rmaxout=300.;
  double *ka=my_malloc(Nk*sizeof(double));
  double *pk=my_malloc(Nk*sizeof(double));
  double *pklin=my_malloc(Nk*sizeof(double));
  double *ra=my_malloc(Nk*sizeof(double));
  double *xi=my_malloc(Nk*sizeof(double));
  double *xilin=my_malloc(Nk*sizeof(double));
  for (int i=0; i<Nk; i++) ka[i]=kmin*pow((kmax/kmin),i*1.0/(Nk-1));
  FILE *fpk, *fxi, *fg;
  char fnamepk[256], fnamexi[256], gbiasfn[256];
  double rsm2=par->r2_smooth+pow(par->l_box/par->n_grid,2)/12.;
  sprintf(gbiasfn,"%s_gbias.txt",par->prefixOut);
  fg=fopen(gbiasfn,"w");
  fprintf (fg,"#1-z 2-r(z) 3-g(z) ");
  for (int ipop=0; ipop<par->n_srcs; ipop++)
    fprintf(fg,"%d-bg_%d(z) ",ipop+4,ipop+1);	
  for (int ipop=0; ipop<par->n_imap; ipop++)
    fprintf(fg,"%d-bi_%d(z) ",ipop+4+par->n_srcs,ipop+1);
  fprintf(fg,"\n");

  // outter loop is over redshifts
  for (double z=0; z<=par->z_max; z+=par->pred_dz) {
    double r=r_of_z(par,z);
    double g=get_bg(par,r,BG_D1,0);
    fprintf (fg,"%g %g %g ",z,r,g);
    for (int i=0; i<Nk; i++) pklin[i]=pk_linear0(par,log10(ka[i]))*g*g;
    pk2xi(Nk,ka,pklin,ra,xilin);
    // inner loop is over populations, ipop=-1 is the unbiased version
#ifdef _DEBUG
    print_info ("Writing predictions of redshift %g:\n",z);
#endif //_DEBUG

    if(par->do_srcs) {
      for (int ipop=0; ipop<par->n_srcs; ipop++) {
	double bias=get_bg(par,r,BG_BZ_SRCS,ipop);
	fprintf(fg,"%g ",bias);
#ifdef _DEBUG
	print_info ("       Population %i, bias %g. \n",ipop,bias);
#endif //_DEBUG
	for (int i=0; i<Nk; i++) pk[i]=pklin[i]*bias*bias*exp(-rsm2*ka[i]*ka[i]);
	pk2xi(Nk,ka,pk,ra,xi);
	for (int i=0; i<Nk; i++) xi[i]=exp(xi[i])-1;
	xi2pk(Nk,ra,xi,ka,pk);
	// now open the files
	sprintf(fnamepk,"%s_pk_srcs_pop%i_z%.3lf.txt",par->prefixOut,ipop,z);
	sprintf(fnamexi,"%s_xi_srcs_pop%i_z%.3lf.txt",par->prefixOut,ipop,z);
	fpk=fopen(fnamepk,"w");
	fprintf (fpk, "# k[h/Mpc] P_tt P_tl P_ll\n");
	fxi=fopen(fnamexi,"w");
	fprintf (fxi, "# r[Mpc/h] xi_tt xi_ll*b^2 xi_ll\n");
	for (int i=0; i<Nk; i++) {
	  if ((ka[i]>=kminout) && (ka[i]<=kmaxout))
	    fprintf (fpk,"%g %g %g %g\n",ka[i],pk[i], pklin[i]*bias*exp(-rsm2*ka[i]*ka[i]), pklin[i]*exp(-rsm2*ka[i]*ka[i]));
	  if ((ra[i]>=rminout) && (ra[i]<=rmaxout))
	    fprintf (fxi,"%g %g %g %g\n",ra[i],xi[i], xilin[i]*bias*bias, xilin[i]);
	}
	fclose(fpk);
	fclose(fxi);
      }
    }
    if(par->do_imap) {
      for (int ipop=0; ipop<par->n_imap; ipop++) {
	double bias=get_bg(par,r,BG_BZ_IMAP,ipop);
	fprintf(fg,"%g ",bias);
#ifdef _DEBUG
	print_info ("       Population %i, bias %g. \n",ipop,bias);
#endif //_DEBUG
	for (int i=0; i<Nk; i++) pk[i]=pklin[i]*bias*bias*exp(-rsm2*ka[i]*ka[i]);
	pk2xi(Nk,ka,pk,ra,xi);
	for (int i=0; i<Nk; i++) xi[i]=exp(xi[i])-1;
	xi2pk(Nk,ra,xi,ka,pk);
	// now open the files
	sprintf(fnamepk,"%s_pk_imap_pop%i_z%.3lf.txt",par->prefixOut,ipop,z);
	sprintf(fnamexi,"%s_xi_imap_pop%i_z%.3lf.txt",par->prefixOut,ipop,z);
	fpk=fopen(fnamepk,"w");
	fprintf (fpk, "# k[h/Mpc] P_tt P_tl P_ll\n");
	fxi=fopen(fnamexi,"w");
	fprintf (fxi, "# k[Mpc/h] xi_tt xi_ll*b^2 xi_ll\n");
	for (int i=0; i<Nk; i++) {
	  if ((ka[i]>=kminout) && (ka[i]<=kmaxout))
	    fprintf (fpk,"%g %g %g %g\n",ka[i],pk[i], pklin[i]*bias, pklin[i]);
	  if ((ra[i]>=rminout) && (ra[i]<=rmaxout))
	    fprintf (fxi,"%g %g %g %g\n",ra[i],xi[i], xilin[i]*bias*bias, xilin[i]);
	}
	fclose(fpk);
	fclose(fxi);
      }
    }
    fprintf(fg,"\n");
  }
  fclose(fg);
  free(ka);
  free(pk);
  free(pklin);
  free(ra);
  free(xi);
  free(xilin);
}
