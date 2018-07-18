#ifndef H__H__SOLVE_SYSTEM__H__H
#define H__H__SOLVE_SYSTEM__H__H

#include "material_properties.h"
#include "hyperelasticity.h"
#include "constitutive_model_handle.h"

#define D_GAMMA_D 0.005
#define D_GAMMA_TOL 1.25

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */
  
extern long perIter_ODE_EXA_metric; //ODE operations accumulated over the current Iter
extern long perTimestep_EXA_metric; //Exa metric accumulated over the current timestep
extern long total_EXA_metric;       //Total Exa metric accumulated over each NR iteration
extern long dof_EXA_metric;         //Exa metric for accumulated Ndof over each NR iteration

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

int staggered_Newton_Rapson_compute(double *pFnp1_out, double *g_out, double *lambda,
                                    double *pFn_in, double *Fn_in, double *Fnp1_in, double *hFn_in, double *hFnp1_in,
                                    double g_n, double dt,
                                    MATERIAL_CONSTITUTIVE_MODEL *mat,
                                    ELASTICITY *elasticity,
                                    CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info,
                                    double *d_gamma,
                                    int *is_it_cnvg);

/// integrate crystal plasticity from pFn to pFnp1 for general cases (includes thermal expansions)
///
/// \param[out] pFnp1_out deformation gradient (plastic) at time step = n + 1
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
int staggered_Newton_Rapson_generalized(double *pFnp1_out, double *g_out, double *lambda, 
                                        double *pFn_in, double *Fn_in, double *Fnp1_in, double *hFn_in, double *hFnp1_in,  
                                        double g_n, double dt, 
                                        MATERIAL_CONSTITUTIVE_MODEL *mat,
                                        ELASTICITY *elasticity, 
                                        CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info);

/// integrate crystal plasticity from pFn to pFnp1
///
/// \param[out] pFnp1_out deformation gradient (plastic) at time step = n + 1
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
int staggered_Newton_Rapson(double *pFnp1_out, double *g_out, double *lambda, 
                            double *pFn_in, double *Fn_in, double *Fnp1_in, 
                            double g_n, double dt, 
                            MATERIAL_CONSTITUTIVE_MODEL *mat,
                            ELASTICITY *elasticity, 
                            CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info); 
                            
int Fnp1_Implicit(double *Fnp1_out, double *Fn_in, double *L_in, double dt);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

class GcmCpIntegrator : public GcmIntegrator
{
  public:

  double gn, gnm1;
  double *gnp1;
  double *lambda;
  double gn_s;
  
  MATERIAL_CONSTITUTIVE_MODEL *mat;
  ELASTICITY *elasticity;
  CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info;

  GcmCpIntegrator(){
    mat         = NULL;
    elasticity  = NULL;
    solver_info = NULL;
  }  
  
  virtual int run_integration_algorithm(const double dt) 
  {
    int err = 0;
    double d_gamma = 0.0;
    int is_it_cnvg = 0;
    
    err += staggered_Newton_Rapson_compute(pFnp1, gnp1, lambda,
                                           pFn_s, Fn_s, F_s, hFn, hFnp1, gn, dt,
                                           mat, elasticity,solver_info,&d_gamma,&is_it_cnvg);
    if((d_gamma)/D_GAMMA_D>D_GAMMA_TOL || !is_it_cnvg)
      return 1;
    
    return err;
                                               
  }
  virtual void set_variable_at_n(void){
    gn_s = gn;
  }
  virtual void update_variable(void){
    gn_s = *gnp1;
  }    
    
};

#endif
