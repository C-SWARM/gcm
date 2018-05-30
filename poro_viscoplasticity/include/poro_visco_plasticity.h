# ifndef H__H__PORO_VISCOPLASTICITY__H__H

#include "material_properties.h"

int poro_visco_plasticity_integration_algorithm(const MaterialPoroViscoPlasticity *mat,
                                                double *Fnp1,
                                                double *Fn,
                                                double *pFnp1,
                                                double *pFn,
                                                double *pc_np1,
                                                const double pc_n,
                                                const double dt); 
# endif
