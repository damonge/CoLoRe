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

void write_predictions(ParamCoLoRe *par) {
  if (!par->do_gals) return;
  if (NodeThis!=0) return;
  print_info("*** Writing predictions (ASCII) \n");
  // first generate k array, sufficiently finely spaced
  int Nk=1000;
  double *ka=my_malloc(Nk*sizeof(double));
  double *pk=my_malloc(Nk*sizeof(double));
  double *ra=my_malloc(Nk*sizeof(double));
  double *xi=my_malloc(Nk*sizeof(double));
  for (int i=0; i<Nk; i++) ka[i]=1e-4*pow((10./1e-4),i*1.0/(Nk-1));
  FILE *fpk, *fxi;
  char fnamepk[256], fnamexi[256];
	
  // outter loop is over redshifts
  for (double z=0; z<=par->z_max; z+=par->pred_dz) {
    // inner loop is over populations, ipop=-1 is the unbiased version
    double r=r_of_z(par,z);
    double g=dgrowth_of_r(par,r);
    for (int ipop=-1; ipop<par->n_gals; ipop++) {
      // first get the linear PS at this redshift
      for (int i=0; i<Nk; i++) pk[i]=pk_linear0(par,log(ka[i]))*g*g;
      // if ipop==-1, just do linear, otherwise do the lognormal prediction
      if (ipop>=0) {
	double bias=bias_of_z_gals(par,z,ipop);
	for (int i=0; i<Nk; i++) pk[i]*=bias*bias*exp(-par->r2_smooth*ka[i]*ka[i]);
	pk2xi(Nk,ka,pk,ra,xi);
	// Fix xi
	for (int i=0; i<Nk; i++) xi[i]=exp(xi[i])-1;
	xi2pk(Nk,ra,xi,ka,pk);
	// now open the files
	sprintf(fnamepk,"%s_pk_pop%i_z%g.dat",par->prefixOut,ipop,z);
	sprintf(fnamexi,"%s_xi_pop%i_z%g.dat",par->prefixOut,ipop,z);

      } else {
	// just calculate xi
	pk2xi(Nk,ka,pk,ra,xi);
	sprintf(fnamepk,"%s_pk_lin_z%g.dat",par->prefixOut,z);
	sprintf(fnamexi,"%s_xi_lin_z%g.dat",par->prefixOut,z);
      }
      printf ("Writing: %s and %s.\n",fnamepk, fnamexi);
      fpk=fopen(fnamepk,"w");
      fxi=fopen(fnamexi,"w");
      for (int i=0; i<Nk; i++) {
	fprintf (fpk,"%g %g \n",ka[i],pk[i]);
	fprintf (fxi,"%g %g \n",ra[i],xi[i]);
      }
      fclose(fpk);
      fclose(fxi);
    }
  }
  free(ka);
  free(pk);
  free(ra);
  free(xi);
}
