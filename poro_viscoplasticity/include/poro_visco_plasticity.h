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
                                       
/// compute conforming pressure for given plastic deformation 
/// \param[in] pJ        determinant of pF
/// \param[in] mat_pvp   poro_viscoplaticity material object
/// \return    computed pc value  
double poro_visco_plasticity_compute_pc(double pJ, 
                                        double pc,
                                        const MaterialPoroViscoPlasticity *mat,
                                        const GcmSolverInfo *solver_info);

void poro_visco_plasticity_update_elasticity_dev(double *eS_out,
                                                 double *L_out,
                                                 const MaterialPoroViscoPlasticity *param,
                                                 double *eF_in,
                                                 const double pc,
                                                 const bool compute_4th_order);

double poro_visco_plasticity_intf_compute_dudj(const double eJ,
                                               const double pc,
                                               const MaterialPoroViscoPlasticity *param);                                               

double poro_visco_plasticity_intf_compute_d2udj2(const double eJ,
                                                 const double pc,
                                                 const MaterialPoroViscoPlasticity *param);                                                 

double poro_visco_plasticity_hardening(const double pc,
                                       const MaterialPoroViscoPlasticity *param);

void poro_visco_plasticity_intf_compute_gammas(double &gamma_dot_d,
                                               double &gamma_dot_v,
                                               double *Fnp1_in,
                                               double *Fn_in,
                                               double *pFnp1_in,
                                               double *pFn_in,
                                               const double pc_np1,
                                               const double pc_n,
                                               const double dt,
                                               const MaterialPoroViscoPlasticity *param,
                                               const GcmSolverInfo *solver_info);
                                             
# endif
