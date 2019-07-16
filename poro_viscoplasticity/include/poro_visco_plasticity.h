/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
/// 
/// [1] University of Notre Dame, Notre Dame, IN

#ifndef H__H__PORO_VISCOPLASTICITY__H__H

#include "material_properties.h"
#include "constitutive_model_handle.h"
#include "GcmSolverInfo.h"

// integration algorithm for PVP model in a staggered manner
int poro_visco_plasticity_integration_algorithm(const MaterialPoroViscoPlasticity *mat,
                                                const GcmSolverInfo *solver_info,
                                                PvpElasticity *elast,
                                                double *Fnp1,
                                                double *Fn,
                                                double *pFnp1,
                                                double *pFn,
                                                double *pc_np1,
                                                const double pc_n,
                                                const double dt_in,
                                                const bool is_implicit = true); 
                                     
/// compute conforming pressure for given plastic deformation 
double poro_visco_plasticity_compute_pc(double pJ, 
                                        double pc,
                                        const MaterialPoroViscoPlasticity *mat,
                                        const GcmSolverInfo *solver_info);

/// compute shear modulus
double poro_visco_plasticity_compute_shear_modulus(const MaterialPoroViscoPlasticity *param,
                                                   const double pc);

/// compute hardening function value of p_c
double poro_visco_plasticity_hardening(const double pc,
                                       const MaterialPoroViscoPlasticity *param);

/// compute plastic strain rates
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

/// compute dMdU for an finite element for PDE
void poro_visco_plasticity_compute_dMdu(double *dMdUs,
                                        double *Grad_us,
                                        const MaterialPoroViscoPlasticity *mat,
                                        const GcmSolverInfo *solver_info,
                                        double *Fnp1_in,
                                        double *Fn_in,
                                        double *pFnp1_in,
                                        double *pFn_in,
                                        const double pc_np1,
                                        const double pc_n,
                                        const double dt,
                                        const int nne,
                                        const int ndofn);

/// compute plastic velocity gradient from PVP model
int poro_visco_plasticity_plastic_velocity_gradient(double *pL_out,
                                                    const MaterialPoroViscoPlasticity *param,
                                                    double *Fnp1,
                                                    double *Fn,
                                                    double *pFnp1,
                                                    double *pFn,
                                                    const double pc_np1,
                                                    const double pc_n);

/// define PVP integrator                                        
class GcmPvpIntegrator : public GcmIntegrator
{
  public:

  double pc_n;
  double *pc_np1;
  double pc_n_s;
  PvpElasticity *elasticity;
  
  MaterialPoroViscoPlasticity *mat;

  GcmPvpIntegrator(){
    pc_np1 = NULL;
    pc_n = pc_n_s = {};
    elasticity = NULL;
  }
  
  GcmPvpIntegrator(MaterialPoroViscoPlasticity *mat_in,
                   PvpElasticity *elasticity_in,
                   GcmSolverInfo *solver_info_in){
    mat         = mat_in;
    elasticity  = elasticity_in;
    solver_info = solver_info_in;
  } 
  
  virtual int constitutive_model_integration(const double dt) 
  {
    int err = 0;
    
    err += poro_visco_plasticity_integration_algorithm(mat, solver_info, elasticity, 
                                                       F_s, Fn_s, pFnp1, pFn_s, 
                                                       pc_np1, pc_n_s, dt);     
    
    return err;
                                               
  }
  virtual void set_variable_at_n(void){
    pc_n_s = pc_n;
  }
  virtual void update_variable(void){
    pc_n_s = *pc_np1;
  }    
    
};                                    

#endif
