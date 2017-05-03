
#include "pvp_interface.h"
#include "KMS-IJSS2017.h"

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

int pvp_intf_perform_integration_algorithm2(double *Fnp1,
                                           double *Fn,
                                           double *pFnp1,
                                           double *pFn,
                                           double *pc_np1,
                                           double pc_n,
                                           KMS_IJSS2017_Parameters *mat_pvp,
                                           double dt)
{
  int err = 0;
/*  
  // Time integrator
  double InitialTime = 0.0; 
  double FinalTime   = dt; 
  double CurrentTime = 0.0;
  size_t TimeStep = 0;
  
  // Parameters and Model definition
  bool Verbose = false;
     
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, NULL, NULL);
    
  ImplicitModelInstance.set_data_from_PDE(Fnp1, Fn, pFnp1, pFn, *pc_np1, pc_n);

  ttl::Tensor<2, DIM_3, double*>  ttl_Fnp1(Fnp1);
  ImplicitModelInstance.StepUpdate(ttl_Fnp1,dt, Verbose);
  ImplicitModelInstance.set_data_to_PDE(pFnp1, pc_np1);
  */
  return err;

}