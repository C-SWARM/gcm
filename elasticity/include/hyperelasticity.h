#ifndef H__H__HYPERELASTICITY__H__H
#define H__H__HYPERELASTICITY__H__H

#include "material_properties.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

struct ELASTICITY;
#ifndef TYPE_ELASTICITY
#define TYPE_ELASTICITY
typedef struct ELASTICITY ELASTICITY;
#endif

typedef int (*elasticity_part)(ELASTICITY *elasticity, double *Fe, int flag);
typedef void (*deviatoric_part)(double *C_in, MATERIAL_ELASTICITY const *mat, double *S);
typedef void (*volume_part)(double *dudj, double J);
typedef int (*compute_elasticity_v1) (ELASTICITY *elasticity, double *V1);
typedef int (*compute_elasticity_v2) (ELASTICITY *elasticity, double *V1, double *V2);

struct ELASTICITY
{
  double *L, *S;
  MATERIAL_ELASTICITY *mat;
  elasticity_part update_elasticity;
  deviatoric_part compute_potential_dev;  
  deviatoric_part compute_PK2_dev;
  deviatoric_part compute_tangent_dev;
  volume_part compute_u;  
  volume_part compute_dudj;
  volume_part compute_d2udj2;
  compute_elasticity_v1 compute_PK2_eff;
  compute_elasticity_v2 compute_Cauchy_eff;
  compute_elasticity_v2 compute_Cauchy;
  
  int compute_stiffness;
};

int construct_elasticity(ELASTICITY *elasticity, MATERIAL_ELASTICITY *mat, const int compute_stiffness);
int destruct_elasticity(ELASTICITY *elasticity);


#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif 
