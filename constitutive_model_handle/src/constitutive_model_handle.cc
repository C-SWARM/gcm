#include <stdio.h>
#include "constitutive_model_handle.h"
#include <ttl/ttl.h>

#define DIM_3   3
#define DIM_3x3 9 

using namespace std;
template<int R=2, int D = 3, class S=double>
using Tensor = ttl::Tensor<R, D, S>;

int constitutive_model_handle_init(CONSTITUTIVE_MODEL_PACKS *cm_pack)
{
  int err = 0;
  cm_pack->mat_e       = NULL;
  cm_pack->mat_p       = NULL;
  cm_pack->mat_d       = NULL;
  cm_pack->solver_info = NULL;
  cm_pack->elasticity  = NULL;
  cm_pack->damage      = NULL;
  cm_pack->cm_mat      = NULL;
  return err;
}

/// find n: n = f(m) where f(
/// L0 = (m-x1)*(m-x2)/(x0-x1)/(x0-x2);
/// L1 = (m-x0)*(m-x2)/(x1-x0)/(x1-x2);
/// L2 = (m-x0)*(m-x2)/(x2-x0)/(x2-x1);
/// 
/// n = y1*L0 + y2*L1 + y3*L2;
/// 
double quadratic_interpolation(const double nm1,
                               const double n,
                               const double np1,
                               const double dt_n,
                               const double dt_np1,
                               const double dt){
                                
  double L0 = dt*(dt - dt_np1)/dt_n/(dt_n+dt_np1);
  double L1 = -(dt + dt_n)*(dt - dt_np1)/dt_n/dt_np1;
  double L2 = (dt + dt_n)*dt/(dt_np1 + dt_n)/dt_np1;

  return nm1*L0 + n*L1 + np1*L2;
}

int GcmIntegrator::run_integration_algorithm(const double dt_n,
                                             const double dt_np1){
  
  int err = 0;
  int stepno = 1;
  
  int not_cnvged = 1;
  
  while(not_cnvged>0)
  {
    double dt_s = dt_np1/stepno;    
        
    // set initial Fn values for sub_cycling
    for(int ib=0; ib<DIM_3x3; ++ib){
      this->Fn_s[ib]  = this->Fn[ib];
      this->pFn_s[ib] = this->pFn[ib];
    }
    this->set_variable_at_n();
          
    // do sub-cycling
    if(solver_info->debug)
      printf("GCM sub-cycling scheme. Number of Sub-cycling = %d\n", stepno);
      
    for(int ia=0; ia<stepno; ++ia)
    {       
      if(solver_info->debug)
        printf(">> Sub-cycling = %d/%d\n", ia, stepno);
      
      // update F(stepno)     
      for(int ib=0; ib<DIM_3x3; ++ib)
      {
        if(stepno==1)
          this->F_s[ib] = this->Fnp1[ib];
        else{
          this->F_s[ib] = quadratic_interpolation(this->Fnm1[ib], 
                                                  this->Fn[ib],
                                                  this->Fnp1[ib],
                                                  dt_n, dt_np1, dt_s*(ia+1));
        }
      }
      
      not_cnvged = this->constitutive_model_integration(dt_s);
      
      if(not_cnvged > 0)
        break;

      for(int ic=0; ic<DIM_3x3; ic++)
      {
        this->pFn_s[ic] = this->pFnp1[ic];
        this->Fn_s[ic] = this->F_s[ic];
      }
      this->update_variable();
    }
    stepno *= 2;
    if(stepno > solver_info->max_subdivision)
      break;
  }
  
  if(not_cnvged>0)
    ++err;
  
  if(err>0 && stepno>1)
    printf("integration algorithm failed after sub-cycling. sub-cycling number %d/%d\n", stepno/2, solver_info->max_subdivision);
  
  return err;
}