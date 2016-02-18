#ifndef H__H__HYPERELASTICITY__H__H
#define H__H__HYPERELASTICITY__H__H

#include "material_properties.h"

struct ELASTICITY;
#ifndef TYPE_ELASTICITY
#define TYPE_ELASTICITY
typedef struct ELASTICITY ELASTICITY;
#endif

typedef int (*elasticity_part)(ELASTICITY *elasticity, double *Fe, int update_stiffness);
typedef void (*deviatoric_part)(double *C_in, MATERIAL_ELASTICITY const *mat, double *S);
typedef void (*volume_part)(double *dudj, double J);

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
  int compute_stiffness;
};

int construct_elasticity(ELASTICITY *elasticity, MATERIAL_ELASTICITY *mat, const int compute_stiffness);
int destruct_elasticity(ELASTICITY *elasticity);

#endif 