#ifndef H__H__SOLVE_SYSTEM__H__H
#define H__H__SOLVE_SYSTEM__H__H

#include "material_properties.h"
#include "hyperelasticity.h"


#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */
  
extern long EXA_metric;

struct CRYSTAL_PLASTICITY_SOLVER_INFO
{
  int max_itr_stag;     // maximum number of iteration for staggerd NR
  int max_itr_hardening; // maximum number of iteration for hardening NR
  int max_itr_M;         // maximum number of iteration for M NR where (M = pFn_I)
  double tol_hardening;
  double tol_M;
  double computer_zero;
  int max_subdivision;
  int debug;
};

typedef struct CRYSTAL_PLASTICITY_SOLVER_INFO CRYSTAL_PLASTICITY_SOLVER_INFO;

int set_crystal_plasticity_solver_info(CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info, 
                                       int max_itr_stag, int max_itr_hardening, int max_itr_M, 
                                       double tol_hardening, double tol_M, double computer_zero);

int print_crystal_plasticity_solver_info(CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info);     

/// integrate crystal plasticity from pFn to pFnp1 for general cases (includes thermal expansions)
///
/// \param[out] pFnp1_out deformation gradient (plastic) at time step = n + 1
/// \param[out] M_out deformation gradient (pFr_I) at time step = n + 1
/// \param[out] g_out updated hardening at time step =  n + 1
/// \param[out] lambda updated Lagrange multiplier at time step =  n + 1
/// \param[in] pFn_in deformation gradient (plastic) at time step = n
/// \param[in] Fn_in deformation gradient (total) at time step = n
/// \param[in] Fnp1_in deformation gradient (total) at time step = n + 1
/// \param[in] hFn_in deformation gradient (due to thermal expansion) at time step = n
/// \param[in] hFnp1_in deformation gradient (due to thermal expansion) at time step = n + 1
/// \param[in] g_n hardening at time step =  n
/// \param[in] dt time step size
/// \param[in] mat material parameters
/// \param[in] elasticity object to compute elasticity (stress)
/// \param[in] solver_info defines numerical parameters 
/// \return non-zero on interal error
int staggered_Newton_Rapson_generalized(double *pFnp1_out, double *M_out, double *g_out, double *lambda, 
                                        double *pFn_in, double *Fn_in, double *Fnp1_in, double *hFn_in, double *hFnp1_in,  
                                        double g_n, double dt, 
                                        MATERIAL_CONSTITUTIVE_MODEL *mat,
                                        ELASTICITY *elasticity, 
                                        CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info);

/// integrate crystal plasticity from pFn to pFnp1
///
/// \param[out] pFnp1_out deformation gradient (plastic) at time step = n + 1
/// \param[out] M_out deformation gradient (pFr_I) at time step = n + 1
/// \param[out] g_out updated hardening at time step =  n + 1
/// \param[out] lambda updated Lagrange multiplier at time step =  n + 1
/// \param[in] pFn_in deformation gradient (plastic) at time step = n
/// \param[in] Fn_in deformation gradient (total) at time step = n
/// \param[in] Fnp1_in deformation gradient (total) at time step = n + 1
/// \param[in] g_n hardening at time step =  n
/// \param[in] dt time step size
/// \param[in] mat material parameters
/// \param[in] elasticity object to compute elasticity (stress)
/// \param[in] solver_info defines numerical parameters 
/// \return non-zero on interal error 
int staggered_Newton_Rapson(double *pFnp1_out, double *M_out, double *g_out, double *lambda, 
                            double *pFn_in, double *Fn_in, double *Fnp1_in, 
                            double g_n, double dt, 
                            MATERIAL_CONSTITUTIVE_MODEL *mat,
                            ELASTICITY *elasticity, 
                            CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info); 
                            
int Fnp1_Implicit(double *Fnp1_out, double *Fn_in, double *L_in, double dt);


#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
