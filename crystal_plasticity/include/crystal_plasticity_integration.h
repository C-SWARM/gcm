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

/// integrate crystal plasticity from pFn to pFnp1
int staggered_Newton_Rapson(double *pFnp1_out, double *g_out, double *lambda,
                            double *pFn_in, double *Fn_in, double *Fnp1_in, double *hFn_in, double *hFnp1_in,
                            double g_n, double dt,
                            MATERIAL_CONSTITUTIVE_MODEL *mat,
                            HyperElasticity *elasticity,
                            GcmSolverInfo *solver_info,
                            double *d_gamma,
                            int *is_it_cnvg);
                            
int Fnp1_Implicit(double *Fnp1_out, double *Fn_in, double *L_in, double dt);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

class GcmCpIntegrator : public GcmIntegrator
{
  public:

  double gn;
  double *gnp1;
  double *lambda;
  double gn_s;
  
  MATERIAL_CONSTITUTIVE_MODEL *mat;
  HyperElasticity *elasticity;

  GcmCpIntegrator(){
    mat         = NULL;
    elasticity  = NULL;
    solver_info = NULL;
    gn = gn_s = {};
    gnp1 = lambda = NULL;
  }
  
  GcmCpIntegrator(MATERIAL_CONSTITUTIVE_MODEL *mat_in,
                  HyperElasticity *elasticity_in,
                  GcmSolverInfo *solver_info_in){

    mat         = mat_in;
    elasticity  = elasticity_in;
    solver_info = solver_info_in;
  } 
  
  virtual int constitutive_model_integration(const double dt)
  {
    int err = 0;
    double d_gamma = 0.0;
    int is_it_cnvg = 0;
    
    err += staggered_Newton_Rapson(pFnp1, gnp1, lambda,
                                   pFn_s, Fn_s, F_s, hFn, hFnp1, gn_s, dt,
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
