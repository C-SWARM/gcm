#include "mpi.h"
#include "pvp_interface.h"
#include "KMS-IJSS2017.h"

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

#define MAX_N_OF_STAGGERED_ITERATIONS 10
#define max_subdivision               16

void print_9x9(double *A,
               FILE *fp,
               const char *name)
{
  fprintf(fp, "double %s[DIM_3x3] = {%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e};\n", name, A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8]);
}
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

int staggered_Newton_Rapson_compute(double *Fnp1,
                                    double *Fn,
                                    double *pFnp1,
                                    double *pFn,
                                    double *pc_np1,
                                    double pc_n,
                                    KMS_IJSS2017_Parameters *mat_pvp,
                                    double dt)
{
  int err = 0;
  
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
  bool Verbose = false;
     
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, NULL, &Time_intg );
    
  ImplicitModelInstance.set_data_from_PDE(Fnp1, Fn, pFnp1, pFn, *pc_np1, pc_n);

  ttl::Tensor<2, DIM_3, double*>  ttl_Fnp1(Fnp1);
  
  try{    
    int itno = ImplicitModelInstance.StepUpdate(ttl_Fnp1,dt, Verbose);
    if(itno >= MAX_N_OF_STAGGERED_ITERATIONS)
      ++err;
  }
  catch(const int itgAlerr)
  {
    ++err;
  }
  ImplicitModelInstance.set_data_to_PDE(pFnp1, pc_np1);
  
  return err;
}

int staggered_Newton_Rapson_subdivision(double *Fnp1,
                                        double *Fn,
                                        double *pFnp1,
                                        double *pFn,
                                        double *pc_np1,
                                        double pc_n,
                                        KMS_IJSS2017_Parameters *mat_pvp,
                                        double dt,
                                        int stepno)
{
  int err = 0;

  ttl::Tensor<2, DIM_3, double> dF,F_s,Fn_s,pFn_s;
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
    
  double dt_s = dt/stepno;

  for(int b=0; b<DIM_3x3; b++)
  {
       dF.data[b] = (Fnp1[b] - Fn[b])/stepno;
     Fn_s.data[b] =  Fn[b];
    pFn_s.data[b] = pFn[b];
  }

  double pc_n_s = pc_n;

  for(int a=0; a<stepno; a++)
  {
    F_s = Fn_s(i,j) + dF(i,j);

    err += staggered_Newton_Rapson_compute(F_s.data, Fn_s.data, pFnp1, pFn_s.data,
                                           pc_np1, pc_n_s, mat_pvp, dt_s);

    if(err > 0)
      break;

    for(int b=0; b<DIM_3x3; b++)
    {
      pFn_s.data[b] = pFnp1[b];
       Fn_s.data[b] = F_s.data[b];
    }
    pc_n_s = *pc_np1;
  }

  return err;
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

  int is_it_cnvg = 0;
  
  is_it_cnvg = staggered_Newton_Rapson_compute(Fnp1, Fn, pFnp1, pFn, pc_np1, pc_n, mat_pvp, dt);

  if(max_subdivision<2)
  {
    if(is_it_cnvg > 0)
      printf("Integration algorithm is diverging\n");
    return is_it_cnvg;
  }

  int stepno = 2;

  while(is_it_cnvg>0)
  {
    is_it_cnvg = staggered_Newton_Rapson_subdivision(Fnp1, Fn, pFnp1, pFn, pc_np1, pc_n, mat_pvp, dt, stepno);
    stepno *= 2;

    if(stepno>max_subdivision)
    {
      ++err;
      break;
    }  
  }

  if(is_it_cnvg>0)
  {
//    int myrank=0;
//    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//    
//    char fname[1024];
//    sprintf(fname, "Fs_%d.txt", myrank); 
//    FILE *fp = fopen(fname, "a");
//        
//    fprintf(fp, "===============================================\n");
//    fprintf(fp, "dt = %e\n",   dt);
//    fprintf(fp, "pc_n = %e\n", pc_n);
//
//    print_9x9( Fnp1, fp, "Fnp1");
//    print_9x9( Fn,   fp, "Fn");
//    print_9x9(pFnp1, fp, "pFnp1");
//    print_9x9(pFn,   fp, "pFn");
//
//    fclose(fp);
      
    printf("After subdivision (%d): Integration algorithm is diverging\n", stepno/2);
  }

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
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp, NULL,NULL);  
    
  return ImplicitModelInstance.compute_dudj(eJ, pc);  
}

double pvp_intf_compute_d2udj2(double eJ,
                               double pc,
                               KMS_IJSS2017_Parameters *mat_pvp)
{
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

void pvp_intf_compute_gammas(double &gamma_d,
                             double &gamma_v,
                             double *Fnp1,
                             double *Fn,
                             double *pFnp1,
                             double *pFn,
                             double pc_np1,
                             double pc_n,
                             KMS_IJSS2017_Parameters *mat_pvp)
{
  
  KMS_IJSS2017_Implicit_BE_Staggered<DIM_3> ImplicitModelInstance(mat_pvp,NULL,NULL);
  ImplicitModelInstance.set_data_from_PDE(Fnp1, Fn, pFnp1, pFn, pc_np1, pc_n);
  
  ImplicitModelInstance.compute_gammas(gamma_d, gamma_v);
}
