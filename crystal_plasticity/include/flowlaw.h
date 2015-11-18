#ifndef H__H__FLOWLAW__H__H
#define H__H__FLOWLAW__H__H

#include "material_properties.h"

double compute_gamma_dot(double *gamma_dots, int N_SYS);

int compute_gamma_dots(double *gamma_dots, double *taus, double g, MATERIAL_CRYSTAL_PLASTICITY *mat_p);

int compute_tau_alphas(double *taus, double *C_in, double *S_in, SLIP_SYSTEM *slip);

int compute_d_gamma_d_tau(double *dgamma_dtaus, double g, double *taus, MATERIAL_CRYSTAL_PLASTICITY *mat_p);

#endif