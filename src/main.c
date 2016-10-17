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

int main(int argc,char **argv)
{ 
  char fnameIn[256];
  ParamCoLoRe *par;
  if(argc!=2) {
    fprintf(stderr,"Usage: ./CoLoRe file_name\n");
    exit(0);
  }
  sprintf(fnameIn,"%s",argv[1]);

  mpi_init(&argc,&argv);

  setbuf(stdout,NULL);
  print_info("\n");
  print_info("|-------------------------------------------------|\n");
  print_info("|                      CoLoRe                     |\n");
  print_info("|-------------------------------------------------|\n\n");

  if(NodeThis==0) timer(4);

  par=read_run_params(fnameIn);

  print_info("Seed : %u\n",par->seed_rng);

  //Create Gaussian density and radial velocity fields
  create_cartesian_fields(par);

  //Interpolate into beams
  pixelize(par);

  //Precompute lensing if needed
  if(par->do_lensing)
    integrate_lensing(par);

  //Poisson-sample the galaxies
  if(par->do_sources)
    get_sources(par);

  //Generate intensity maps
  if(par->do_imap)
    get_imap(par);

  //Generate kappa maps
  if(par->do_kappa)
    get_kappa(par);

  //Write output
  if(par->do_sources)
    write_catalog(par);
  if(par->do_imap)
    write_imap(par);
  if(par->do_kappa)
    write_kappa(par);
  if(par->do_pred)
    write_predictions(par);

  print_info("\n");
  print_info("|-------------------------------------------------|\n\n");

  param_colore_free(par);

  if(NodeThis==0) timer(5);
#ifdef _HAVE_MPI
  MPI_Finalize();
#endif //_HAVE_MPI
  return 0;
}
