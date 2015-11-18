#ifndef H__H__STRAIN_ENERGY_DENSITY_FUNCTION__H__H
#define H__H__STRAIN_ENERGY_DENSITY_FUNCTION__H__H

#include "material_properties.h"

void SEDF_devStress_Mooney_Rivlin(double *C_in,
			     MATERIAL_ELASTICITY const *mat,
			     double *S);

void SEDF_devStress_Linear(double *C,
		      MATERIAL_ELASTICITY const *mat,
		      double *S);
		      
void SEDF_matStiff_Mooney_Rivlin(double *C_in,
			    MATERIAL_ELASTICITY const *mat,
			    double *L_out);
			    
void SEDF_matStiff_Linear(double *C,
		     MATERIAL_ELASTICITY const *mat,
		     double *L);

void SEDF_dUdJ_Common(double *dUdJ, double J);

void SEDF_dUdJ_Doll_Schweizerhof_7(double *dUdJ, double J);

void SEDF_dUdJ_Doll_Schweizerhof_8(double *dUdJ, double J);

void SEDF_d2UdJ2_Common_new(double *d2UdJ2, double const J);

void SEDF_d2UdJ2_Doll_Schweizerhof_7(double *d2UdJ2, double J);

void SEDF_d2UdJ2_Doll_Schweizerhof_8(double *d2UdJ2, double J);

#endif
