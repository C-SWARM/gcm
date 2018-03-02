
#include "pvp_interface.h"
#include "KMS-IJSS2017.h"

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

double pvp_intf_hardening_law(double pc,
                              KMS_IJSS2017_Parameters *mat_pvp)
{
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, NULL,NULL);
  return ImplicitModelInstance.HardeningLaw(pc);
}

/// compute conforming pressure for given plastic deformation 
/// \param[in] pJ        determinant of pF
/// \param[in] mat_pvp   poro_viscoplaticity material object
/// \return    computed pc value  
double pvp_intf_compute_pc(double pJ, double pc,
                           KMS_IJSS2017_Parameters *mat_pvp)
{
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, NULL,NULL);
  double logJp = log(pJ);
  ImplicitModelInstance.FindpcFromJp(logJp, pc, 1.0e-12, false );
        
  return pc;
}

int pvp_intf_perform_integration_algorithm(double *Fnp1,
                                           double *Fn,
                                           double *pFnp1,
                                           double *pFn,
                                           double *pc_np1,
                                           double pc_n,
                                           KMS_IJSS2017_Parameters *mat_pvp,
                                           double dt)
{
  int err = 0;
  unsigned PrintEveryNSteps= 100;

  std::string pathstr = "./out/KMS_IJSS2017";
  std::string logstr = pathstr + ".implicit.cst.smb.log";
  std::string Fstr = pathstr + ".implicit.cst.smb.F.txt";
  std::string pFstr = pathstr + ".implicit.cst.smb.Fp.txt";
  std::string Sstr = pathstr + ".implicit.cst.smb.S.txt";
  std::string KSstr = pathstr + ".implicit.cst.smb.KS.txt";
  std::string sigmastr = pathstr + ".implicit.cst.smb.sigma.txt";
  std::string pcstr = pathstr + ".implicit.cst.smb.pc.txt";
  std::string matstr = pathstr + ".implicit.cst.smb.mat.txt";
    
//  KMS_IJSS2017_IO IO(logstr, Fstr, pFstr, Sstr, KSstr, sigmastr, pcstr, matstr, PrintEveryNSteps);
    
  
  // Time integrator
  double InitialTime = 0.0; 
  double FinalTime   = dt; 
  double CurrentTime = 0.0;
  size_t TimeStep = 0;
  
  TimeIntegrationDataManager Time_intg(InitialTime, FinalTime, dt, CurrentTime, TimeStep);
  
  

  // Parameters and Model definition
  bool usingSmoothMacauleyBrackets = true;
  bool Verbose = false;
     
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, NULL, &Time_intg );
    
  ImplicitModelInstance.set_data_from_PDE(Fnp1, Fn, pFnp1, pFn, *pc_np1, pc_n);

  ttl::Tensor<2, DIM_3, double*>  ttl_Fnp1(Fnp1);
  
  try{
    ImplicitModelInstance.StepUpdate(ttl_Fnp1,dt, Verbose);    
  }
  catch(const int itgAlerr)
  {
    ++err;
  }
  ImplicitModelInstance.set_data_to_PDE(pFnp1, pc_np1);
  
  return err;
}

int pvp_intf_update_elasticity(double *eF,
                               double pc,
                               double *S_in,
                               double *L_in,
                               KMS_IJSS2017_Parameters *mat_pvp,
                               int compute_stiffness)
{
  int err = 0;
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, NULL,NULL);  
    
  ImplicitModelInstance.update_elasticity(eF, pc, S_in, L_in, compute_stiffness);  

  return err;
}

int pvp_intf_update_elasticity_dev(double *eF,
                                   double pc,
                                   double *S_in,
                                   double *L_in,
                                   KMS_IJSS2017_Parameters *mat_pvp,
                                   int compute_stiffness)
{
  int err = 0;
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, NULL,NULL);  
    
  ImplicitModelInstance.update_elasticity_dev(eF, pc, S_in, L_in, compute_stiffness);  

  return err;
}

double pvp_intf_compute_dudj(double eJ,
                             double pc,
                             KMS_IJSS2017_Parameters *mat_pvp)
{
  int err = 0;
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, NULL,NULL);  
    
  return ImplicitModelInstance.compute_dudj(eJ, pc);  
}

double pvp_intf_compute_d2udj2(double eJ,
                               double pc,
                               KMS_IJSS2017_Parameters *mat_pvp)
{
  int err = 0;
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, NULL,NULL);  
    
  return ImplicitModelInstance.compute_d2udj2(eJ, pc);  
}

int pvp_intf_compute_dMdF(double *dMdF_in,
                          double *Fnp1,
                          double *Fn,
                          double *pFnp1,
                          double *pFn,
                          double pc_np1,
                          double pc_n,
                          KMS_IJSS2017_Parameters *mat_pvp,
                          double dt)
                          
{
  int err = 0;
  
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp,NULL,NULL);
  ImplicitModelInstance.set_data_from_PDE(Fnp1, Fn, pFnp1, pFn, pc_np1, pc_n);
  
  ttl::Tensor<4, DIM_3, double*> dMdF(dMdF_in);
  ImplicitModelInstance.compute_dMdF(dMdF.data, dt);
  
  return err;
}
