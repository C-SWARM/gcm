/**
 * Authors:
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */
 
# ifndef H__H__PORO_VISCOPLASTICITY__H__H

#include "material_properties.h"
#include "GcmSolverInfo.h"

int poro_visco_plasticity_integration_algorithm(const MaterialPoroViscoPlasticity *mat,
                                                const GcmSolverInfo *solver_info,
                                                double *Fnp1,
                                                double *Fn,
                                                double *pFnp1,
                                                double *pFn,
                                                double *pc_np1,
                                                const double pc_n,
                                                const bool is_implicit = true); 

void poro_visco_plasticity_update_elasticity(double *eS_out,
                                             double *L_out,
                                             const MaterialPoroViscoPlasticity *param,
                                             double *eF_in,
                                             const double pc,
                                             const bool compute_4th_order = false);

double poro_visco_plasticity_hardening(const double pc,
                                       const MaterialPoroViscoPlasticity *param);
                                             
# endif
