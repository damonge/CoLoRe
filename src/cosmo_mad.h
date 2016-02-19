///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CosmoMAD.                                    //
//                                                                   //
// CosmoMAD is free software: you can redistribute it and/or modify  //
// it under the terms of the GNU General Public License as published //
// by the Free Software Foundation, either version 3 of the License, //
// or (at your option) any later version.                            //
//                                                                   //
// CosmoMAD is distributed in the hope that it will be useful, but   //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CosmoMAD.  If not, see <http://www.gnu.org/licenses/>. //
//                                                                   //
///////////////////////////////////////////////////////////////////////
#ifndef _COSMO_MAD_H
#define _COSMO_MAD_H

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>

#define CSM_FOURPITHIRD 4.18879020478639 //4*pi/3
#define CSM_LOGTEN 2.302585092994 //ln(10)
#define CSM_RTOD 57.29577951308232 //180/pi
#define CSM_DTOR 0.017453292519943295 //pi/180
#define CSM_HGYR 9.777751486751187 //H0^-1 in Gyrs/h
#define CSM_HMPC 2997.92458 //H0^-1 in Mpc/h

#define CSM_MIN(a, b)  (((a) < (b)) ? (a) : (b))

typedef struct {
  //Tuneable parameters
  double OM;
  double OB;
  double OL;
  double h;
  double w0;
  double wa;
  double TCMB;

  // Derived parameters
  double OK;
  int ksign;
  int normalDE;
  int constantw;
  double phorizon;
  double a_equality;
  double growth0;

  // a(t) spline
  int at_spline_set;
  gsl_interp_accel *intacc_at;
  gsl_spline *spline_at;
} Csm_bg_params;

typedef struct {
  /* Background parameters */
  int bg_params_set;
  Csm_bg_params *bg;
} Csm_params;


void csm_unset_gsl_eh(void);
void csm_set_verbosity(int verb);
void csm_params_free(Csm_params *par);
Csm_params *csm_params_new(void);
double csm_hubble(Csm_params *par,double aa);
double csm_particle_horizon(Csm_params *par,double aa);
double csm_radial_comoving_distance(Csm_params *par,double aa);
double csm_curvature_comoving_distance(Csm_params *par,double aa);
double csm_growth_factor(Csm_params *par,double aa);
double csm_f_growth(Csm_params *par,double aa);
void csm_background_set(Csm_params *par,
			double OmegaM,double OmegaL,double OmegaB,
			double ww,double wwa,double hh,double T_CMB);

#endif //_COSMO_MAD_
