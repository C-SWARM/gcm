#ifndef H__H__HARDENING__H__H
#define H__H__HARDENING__H__H

#include "material_properties.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

double compute_g_np1(double gs_np1, double gn, double dt, double gamma_dot_np1, 
                     MATERIAL_CRYSTAL_PLASTICITY *mat_p);

double compute_gs_np1(double gamma_dot_np1, MATERIAL_CRYSTAL_PLASTICITY *mat_p);


#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
