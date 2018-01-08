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
  int test_memory=0;
  char fnameIn[256];
  ParamCoLoRe *par;
  if(argc<2) {
    fprintf(stderr,"Usage: ./CoLoRe file_name\n");
    exit(0);
  }
  if(!strcmp(argv[1],"--test-memory")) {
    test_memory=1;
    sprintf(fnameIn,"%s",argv[2]);
  }
  else
    sprintf(fnameIn,"%s",argv[1]);

  mpi_init(&argc,&argv);

  setbuf(stdout,NULL);
  print_info("\n");
  print_info("|-------------------------------------------------|\n");
  print_info("|                      CoLoRe                     |\n");
  print_info("|-------------------------------------------------|\n\n");

  if(NodeThis==0) timer(4);

  par=read_run_params(fnameIn,test_memory);
  if(par==NULL) {
    if(NodeThis==0) timer(5);
#ifdef _HAVE_MPI
    MPI_Finalize();
#endif //_HAVE_MPI
    return 0;
  }

  print_info("Seed : %u\n",par->seed_rng);

  if(par->do_pred)
    write_predictions(par);
  if(!(par->just_do_pred)) {

    //Create Gaussian density and radial velocity fields
    create_cartesian_fields(par);
    
    //Lognormalize density field
    compute_physical_density_field(par);
    
    //Compute normalization of density field for biasing
    compute_density_normalization(par);
    
    //Get information from slabs
    if(par->do_srcs)
      srcs_set_cartesian(par);
    if(par->do_imap)
      imap_set_cartesian(par);
    if(par->do_kappa)
      kappa_set_cartesian(par);
    if(par->do_isw)
      isw_set_cartesian(par);
    
    //Distribute information across
    if(par->do_srcs)
      srcs_distribute(par);
    if(par->do_imap)
      imap_distribute(par);
    if(par->do_kappa)
      kappa_distribute(par);
    if(par->do_isw)
      isw_distribute(par);
    
    //Postprocess after 
    if(par->do_srcs)
      srcs_get_local_properties(par);
    if(par->do_imap)
      imap_get_local_properties(par);
    if(par->do_kappa)
      kappa_get_local_properties(par);
    if(par->do_isw)
      isw_get_local_properties(par);
    
    //All-to-all communication of density field
    //and computation of all required quantities
    if(par->need_beaming)
      get_beam_properties(par);
    
    //Write output
    if(par->do_srcs)
      write_srcs(par);
    if(par->do_imap)
      write_imap(par);
    if(par->do_kappa)
      write_kappa(par);
    if(par->do_isw)
      write_isw(par);
  }

  print_info("\n");
  print_info("|-------------------------------------------------|\n\n");
  
  param_colore_free(par);

  if(NodeThis==0) timer(5);
#ifdef _HAVE_MPI
  MPI_Finalize();
#endif //_HAVE_MPI
  return 0;
}
