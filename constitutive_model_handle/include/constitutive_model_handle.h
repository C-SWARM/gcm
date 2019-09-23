#ifndef H__H__CONSTITUTIVE_MODEL_HANDLE__H__H
#define H__H__CONSTITUTIVE_MODEL_HANDLE__H__H

#include<stdlib.h>
#include <vector>
#include"GcmSolverInfo.h"
#include "hyperelasticity.h"

extern long perIter_ODE_EXA_metric; //ODE operations accumulated over the current Iter
extern std::vector<long> ODE_EXA_metric;         //Exa metric accumulated over the current timestep
extern std::vector<long> dof_EXA_metric;         //Exa metric for accumulated Ndof over each NR iteration

struct MATERIAL_ELASTICITY;
#ifndef TYPE_MATERIAL_ELASTICITY
#define TYPE_MATERIAL_ELASTICITY
typedef struct MATERIAL_ELASTICITY MATERIAL_ELASTICITY;
#endif

struct MATERIAL_CRYSTAL_PLASTICITY;
#ifndef TYPE_MATERIAL_CRYSTAL_PLASTICITY
#define TYPE_MATERIAL_CRYSTAL_PLASTICITY
typedef struct MATERIAL_CRYSTAL_PLASTICITY MATERIAL_CRYSTAL_PLASTICITY;
#endif

struct MATERIAL_CONTINUUM_DAMAGE;
#ifndef TYPE_MATERIAL_CONTINUUM_DAMAGEY
#define TYPE_MATERIAL_CONTINUUM_DAMAGE
typedef struct MATERIAL_CONTINUUM_DAMAGE MATERIAL_CONTINUUM_DAMAGE;
#endif

struct MATERIAL_J2_PLASTICITY;
#ifndef TYPE_MATERIAL_J2_PLASTICITY
#define TYPE_MATERIAL_J2_PLASTICITY
typedef struct MATERIAL_J2_PLASTICITY MATERIAL_J2_PLASTICITY;
#endif

struct CRYSTAL_PLASTICITY_SOLVER_INFO;
#ifndef TYPE_CRYSTAL_PLASTICITY_SOLVER_INFO
#define TYPE_CRYSTAL_PLASTICITY_SOLVER_INFO
typedef struct CRYSTAL_PLASTICITY_SOLVER_INFO CRYSTAL_PLASTICITY_SOLVER_INFO;
#endif

struct CONTINUUM_DAMAGE;
#ifndef TYPE_CONTINUUM_DAMAGE
#define TYPE_CONTINUUM_DAMAGE
typedef struct CONTINUUM_DAMAGE CONTINUUM_DAMAGE;
#endif

struct MATERIAL_CONSTITUTIVE_MODEL;
#ifndef TYPE_MATERIAL_CONSTITUTIVE_MODEL
#define TYPE_MATERIAL_CONSTITUTIVE_MODEL
typedef struct MATERIAL_CONSTITUTIVE_MODEL MATERIAL_CONSTITUTIVE_MODEL;
#endif

struct CONSTITUTIVE_MODEL_PACKS;
#ifndef TYPE_CONSTITUTIVE_MODEL_PACKS
#define TYPE_CONSTITUTIVE_MODEL_PACKS
typedef struct CONSTITUTIVE_MODEL_PACKS CONSTITUTIVE_MODEL_PACKS;
#endif

struct CONSTITUTIVE_MODEL_PACKS 
{
  MATERIAL_ELASTICITY            *mat_e;
  MATERIAL_CRYSTAL_PLASTICITY    *mat_p;
  MATERIAL_CONTINUUM_DAMAGE      *mat_d;
  MATERIAL_J2_PLASTICITY         *mat_j2p;
  CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info;
  HyperElasticity                *elasticity;
  CONTINUUM_DAMAGE               *damage;
  MATERIAL_CONSTITUTIVE_MODEL    *cm_mat;        
};


int constitutive_model_handle_init(CONSTITUTIVE_MODEL_PACKS *cm_pack);

/// define generalized constitutive model integrator
class GcmIntegrator
{
  public:
    GcmSolverInfo *solver_info;
    
    double *Fnp1, *Fn, *Fnm1, *pFnp1, *pFn;
    double *hFnp1, *hFn;
    double F_s[9], Fn_s[9], pFn_s[9];
    GcmIntegrator(){ 
      Fnp1 = Fn = Fnm1 = pFnp1 = pFn = NULL;
      for(int ia=0; ia<9; ia++)
        F_s[ia] = Fn_s[ia] = pFn_s[ia] = 0.0;
    }

    int run_integration_algorithm(const double dt_n,
                                  const double dt_np1);
    
    virtual int constitutive_model_integration(const double dt){ return 0; };

    virtual void set_variable_at_n(void){};
    
    virtual void update_variable(void){};
    
    void set_tensors(double *Fnp1_in,
                     double *Fn_in,
                     double *Fnm1_in,
                     double *pFnp1_in,
                     double *pFn_in,
                     double *hFnp1_in,
                     double *hFn_in){
      Fnp1  = Fnp1_in;
      Fn    = Fn_in;
      Fnm1  = Fnm1_in;
      pFnp1 = pFnp1_in;
      pFn   = pFn_in;
      hFnp1 = hFnp1_in;
      hFn   = hFn_in;
    }                   
};


/// find n: n = f(m) where f(m) is given below:
/// L0 = (m-x1)*(m-x2)/(x0-x1)/(x0-x2);
/// L1 = (m-x0)*(m-x2)/(x1-x0)/(x1-x2);
/// L2 = (m-x0)*(m-x2)/(x2-x0)/(x2-x1);
/// 
/// n = y1*L0 + y2*L1 + y3*L2;
/// 
double quadratic_interpolation(const double nm1,
                               const double n,
                               const double np1,
                               const double dtn,
                               const double dtnp1,
                               const double dt);
                              
#endif
