/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
/// 
/// [1] University of Notre Dame, Notre Dame, IN

#include <iostream>
#include <iomanip>
#include <math.h>
#include <sys/time.h> // clock_t, clock, CLOCKS_PER_SEC
#include <string.h>
#include <ttl/ttl.h>

#include "input_file_read_handles.h"
#include "hyperelasticity.h"
#include "crystal_plasticity_integration.h"
#include "GcmSolverInfo.h"

#include"poro_visco_plasticity.h"

using namespace std;
template<int R=2, int D = 3, class S=double>
using Tensor = ttl::Tensor<R, D, S>;

static constexpr ttl::Index<'i'> i;
static constexpr ttl::Index<'j'> j;
static constexpr ttl::Index<'k'> k;
static constexpr ttl::Index<'l'> l;
    
#define DIM_3x3 9

class SIM_PARAMS
{
public:
  double dt;             /// time step size
  int    stepno;         /// number of time steps
  char   file_out[1024]; /// output file name
  int    intg_type;      ///
  int    Lno;            /// number velocity gradient
  double velocity;
  int    dim; 
  double *t_end;
  double *L;
  bool   implicit;
  double TMD;
  SIM_PARAMS()
  {
    dt        = 0.0;
    stepno    = 0;
    intg_type = 0;
    Lno       = 0;
    velocity  = 0.0;
    dim       = 0;
    t_end = NULL;
    L     = NULL;
  }
  ~SIM_PARAMS()
  {
    if(t_end!=NULL)
      delete[] t_end;
      
    if(L !=NULL)
      delete[] L; 
  }
};


enum param_names {
  PARAM_yf_M,       // Yield function parameters
  PARAM_yf_alpha,   //   :
  PARAM_flr_m,      // Flow rule parameters
  PARAM_flr_gamma0, //   :
  PARAM_hr_a1,      // Hardening rule parameters
  PARAM_hr_a2,      //   :
  PARAM_hr_Lambda1, //   :
  PARAM_hr_Lambda2, //   :
  PARAM_c_inf,      // Cohesion rule parameters
  PARAM_c_Gamma,    //   :
  PARAM_d_B,        // Transition rule parameters
  PARAM_d_pcb,      //   :
  PARAM_mu_0,       // Shear modulus parameters
  PARAM_mu_1,       //   :
  PARAM_K_p0,       // Bulk modulus parameters
  PARAM_K_kappa,    //   :
  PARAM_pl_n,       // Power law exponent
  PARAM_cf_g0,      // Compaction function parameters
  PARAM_cf_pcinf,   //   :
  PARAM_TMD,        // initial TMD, TDM*pJ = TMD_0;
  PARAM_NO
};

int print_results(FILE *fp,
                  double *F_in,
                  double *pF_in,
                  const double pc,
                  const double t,
                  const MaterialPoroViscoPlasticity *mat_pvp,
                  const double *L0,
                  const int print_option){
  int err = 0;
  
  Tensor<2,3,double*> F(F_in), pF(pF_in);
  Tensor<2,3,double> eF, pF_I;
  pF_I = ttl::inverse(pF);
  eF = F(i,k)*pF_I(k,j);
  double eJ = ttl::det(eF);
  double pJ = ttl::det(pF);
  double  J = ttl::det(F);

  Tensor<2,3,double> S, sigma;
  poro_visco_plasticity_update_elasticity(S.data, NULL, mat_pvp, eF.data, pc, false);
  sigma(i,l) = eF(i,j)*S(j,k)*eF(l,k)/eJ;

  if(print_option==1)        
    printf("%.17e %.17e %.17e\n", t, pc, pJ);

  fprintf(fp, "%.17e %.17e %.17e %.17e %.17e ", t, pc, J, eJ, pJ);

  for(int ia=0; ia<DIM_3x3; ia++)
    fprintf(fp, "%.17e ", sigma.get(ia));

  double l[3]={};
  for(int ia=0; ia<3; ia++){
    l[ia] = F_in[ia*3+0] + F_in[ia*3+1] + F_in[ia*3+2];
    fprintf(fp, "%.17e ", (l[ia]-L0[ia])/L0[ia]);
  }
      
  fprintf(fp, "\n");
  
  return err;
}
                  
/// compute total deformation gradient by integrating velocity gradient
///
/// \param[out] *F   deformation gradient at t
/// \param[in]   t   current time
/// \return non-zero on internal error
int F_of_t(double *F,
           double t);
            

/// compute total deformation gradient by integrating velocity gradient
///
/// \param[in]  *Fn  deformation gradient at tn
/// \param[out] *F   deformation gradient at t
/// \param[out] *L   velocity gradient
/// \param[in]   t   current time
/// \param[in]  &sim simulation parameters are defined in this object e.g. dt and time step size and velocity gradient 
/// \return non-zero on internal error
int F_of_t(double *Fn,
           double *F,
           double *L,
           double t,            
           const SIM_PARAMS &sim)
{
  int err = 0;

  switch(sim.intg_type)
  {
    case 0:
    {
      int ia = 0;
      for(; ia<sim.Lno; ia++)
      {
        if(t<sim.t_end[ia])
          break;
      }
      memcpy(L,sim.L+ia,DIM_3x3*sizeof(double));
      Fnp1_Implicit(F, Fn, L, sim.dt);      
      break;
    }
    case 1:
    {
      double I[DIM_3x3] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};      
      memcpy(F, I, DIM_3x3*sizeof(double));
      if(sim.dim>=1) F[0] = 1.0 + sim.velocity*t;
      if(sim.dim>=2) F[4] = 1.0 + sim.velocity*t;
      if(sim.dim>=3) F[8] = 1.0 + sim.velocity*t;        
      break;      
    }
    case 2:
      err += F_of_t(F, t);
      break;
       
    default:
      break;
  }      
  

  return err;
};

/// read input file
///
/// Input file is formated as:
/// # is comments
/// 1 # if 1: implicit 
/// #      0: explicit
/// #-------------------------------------------------------------------------
/// #  yf_M  yf_alpha  flr_m flr_gamma0  hr_a1 hr_a2 hr_Lambda1 hr_Lambda2  c_inf c_Gamma d_B d_pcb mu_0 mu_1 K_p0  K_kappa pl_n cf_g0 cf_pcinf TDM
/// #-------------------------------------------------------------------------
///    1.0   1.1       0.15  0.0005      0.62  0.37  77.22      13.01       15    0.01    0.2 5.8   30   60   0.063 0.008   2    1     290 0.9
/// #-------------------------------------------------------------------------
/// # Analysis name
/// #-------------------------------------------------------------------------
///   iso_compaction
/// #
/// #-------------------------------------------------------------------------
/// # time steps
/// # dt     number_of_time_steps
/// #-------------------------------------------------------------------------
///   0.01 13200
/// #-------------------------------------------------------------------------
/// # integration type (method of computing deformation gradient)
/// # 0: F_of_t(L1,L2,L3, ...) : L = n number of velocity gradient with t(end) to be transient when t reaches t_end
/// # 1: F_of_t(v,dim)         : v = velocity (F = 1 + v*t), dim = 1: uniaxial displacement
/// #                                                              2: biaxial displacement
/// #                                                              3: triaxial displacement
/// # 2: user define F_of_t    : requres compile
/// #-------------------------------------------------------------------------
/// # [integration type] if type==0: number of Ls and list of Ls are followed]
/// #                    if type==1: velocity and dim are followed
/// 0 3
/// # 1 50e-6 1
/// #-------------------------------------------------------------------------
/// # 1st velocity gradient
/// # L11 L12 L13 L21 L22  L23 L31 L32 L33  t_end
/// #-------------------------------------------------------------------------
///   -0.005  0.0    0.0
///    0.0   -0.005  0.0
///    0.0    0.0   -0.005 42.6
/// 
///   -0.005  0.005  0.0
///    0.0   -0.005  0.0
///    0.0    0.0   -0.005 46.0
/// 
///    0.0    0.005  0.0
///    0.0    0.0    0.0
///    0.0    0.0    0.0 132.0 
/// \param[in]      *in      file pointer to be read
/// \param[in, out] &mat_pvp material property object for porovisco-plasticity
/// \param[in, out] &sim     simulation parameter object such as numer of time step, time step size
/// \return         non-zero on interal error
int read_input_file(const char *input_file,
                    MaterialPoroViscoPlasticity &mat_pvp,
                    SIM_PARAMS &sim)
{
  int err = 0;
  
  FILE *fp_sim = NULL;
  fp_sim = fopen(input_file, "r");

  // read material and simulation parameters
  // if fail to read the file, exit    
  if(fp_sim==NULL)
  { 
    cout << "fail to read ["<< input_file << "]\n";
    exit(-1);      
  }  

  err += goto_valid_line(fp_sim);
  
  int implicit;
  fscanf(fp_sim, "%d", &implicit);
  
  sim.implicit = true;
  if(implicit==0)
    sim.implicit = false;
    
  err += goto_valid_line(fp_sim);
      
  // read material parameters 
  double param[PARAM_NO];
  for(int ia=0; ia<PARAM_NO; ia++)
    fscanf(fp_sim, "%lf", param+ia);
    
  set_properties_poro_visco_plasticity(&mat_pvp,
                                       param[PARAM_yf_M],
                                       param[PARAM_yf_alpha],
                                       param[PARAM_flr_m],
                                       param[PARAM_flr_gamma0],
                                       param[PARAM_hr_a1],
                                       param[PARAM_hr_a2],
                                       param[PARAM_hr_Lambda1],
                                       param[PARAM_hr_Lambda2],
                                       param[PARAM_c_inf],
                                       param[PARAM_c_Gamma],
                                       param[PARAM_d_B],
                                       param[PARAM_d_pcb],
                                       param[PARAM_mu_0],
                                       param[PARAM_mu_1],
                                       param[PARAM_K_p0],
                                       param[PARAM_K_kappa],
                                       param[PARAM_pl_n],
                                       param[PARAM_cf_g0],
                                       param[PARAM_cf_pcinf]);

  sim.TMD = param[PARAM_TMD];
  
  err += goto_valid_line(fp_sim);
  fscanf(fp_sim, "%s", sim.file_out);
  
  err += goto_valid_line(fp_sim);
  fscanf(fp_sim, "%lf %d", &(sim.dt), &(sim.stepno));
  
  err += goto_valid_line(fp_sim);
  fscanf(fp_sim, "%d", &(sim.intg_type));
  
  switch(sim.intg_type)
  {
    case 0:
    {  
      fscanf(fp_sim, "%d", &(sim.Lno));
      sim.L = new double[sim.Lno*DIM_3x3];
      sim.t_end = new double[sim.Lno];
      
      // read 1st velocity gradient
      for(int ia=0; ia<sim.Lno; ia++)
      {
        err += goto_valid_line(fp_sim);
        double *L = sim.L + DIM_3x3*ia;
        fscanf(fp_sim, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", L+0,L+1,L+2,
                                                                  L+3,L+4,L+5,
                                                                  L+6,L+7,L+8,
                                                                  sim.t_end + ia);
      }
      break;
    }
    case 1:
    {
      fscanf(fp_sim, "%lf %d", &(sim.velocity), &(sim.dim));
      break;
    }
    case 2:
      // use user define F_of_t
      // do_nothing
      break;
    default:
      // use user define F_of_t
      // do_nothing
      break;
  }
  fclose(fp_sim);
  
  return err;
}                    

/// print numerical and material parameters
/// 
/// \param[in, out] &mat_pvp material property object for porovisco-plasticity
/// \param[in, out] &sim     simulation parameter object such as numer of time step, time step size
/// \return         non-zero on interal error 
int print_inputs(MaterialPoroViscoPlasticity &mat_pvp,
                 SIM_PARAMS &sim)
{
  int err = 0;

  cout << "--------------------------------------------" << endl;
  cout << "Numerical parameters" << endl;  
  cout << "--------------------------------------------" << endl;
  if(sim.implicit)
    cout << "integration scheme\t: implicit" << endl;
  else
    cout << "integration scheme\t: explicit" << endl;
            
  cout << "output file name\t: "     << sim.file_out  << endl;
  cout << "dt\t\t\t: "               << sim.dt        << endl;
  cout << "number of time steps\t: " << sim.stepno    << endl;
  cout << "integration type\t: "     << sim.intg_type << endl;
  cout << "loading velocity \t: "    << sim.velocity  << endl;
  cout << "loading dimension \t: "   << sim.dim       << endl;  
  cout << "number of Ls\t\t: "       << sim.Lno       << endl;
  cout << "TMD         \t\t: "       << sim.TMD       << endl;
  for(int ia=0; ia<sim.Lno; ia++)
  {
    cout << "End time = " << sim.t_end[ia] << "\t: ";
    cout << "L" << "(" << ia << ") = ["; 
    cout << sim.L[ia*DIM_3x3+0] << " " << sim.L[ia*DIM_3x3+1] << " " << sim.L[ia*DIM_3x3+2] << "; ";
    cout << sim.L[ia*DIM_3x3+3] << " " << sim.L[ia*DIM_3x3+4] << " " << sim.L[ia*DIM_3x3+5] << "; ";
    cout << sim.L[ia*DIM_3x3+6] << " " << sim.L[ia*DIM_3x3+7] << " " << sim.L[ia*DIM_3x3+8] << "];" << endl;
  }
  
  cout << "--------------------------------------------" << endl;
  cout << "Material parameters" << endl;  
  cout << "--------------------------------------------" << endl;
  print_material_property_poro_visco_plasticity(&mat_pvp);

  return err;
}                 

int main(int argc,char *argv[])
{
  int err = 0;
    
  char fn_sim[1024];
  int print_option = 0;

  if(argc<2)
  {
    cout << "Usage: mpirun -np # ./poro_viscoplasticity [FILE_SIM] [print_option]\n";
    cout << "\t[FILE_SIM]\t: file path with simulation parameters\n";
    cout << "\t[print_option]\t: if -1: do not print anything\n";
    cout << "\t                  if  0: print input parameters\n";
    cout << "\t                  if  1: print input and outputs\n";
    cout << "\t                  default is 0\n";    
    exit(-1);
  }
  else
  {    
    // read commend line arguments
    sscanf(argv[1], "%s", fn_sim);
    
    if(argc>=3)
      sscanf(argv[2], "%d", &print_option);
    
    MaterialPoroViscoPlasticity mat_pvp;
    SIM_PARAMS sim;

    err += read_input_file(fn_sim, mat_pvp, sim);
    
    GcmSolverInfo solver_info;
    set_gcm_solver_info(&solver_info, 10, 10, 100, 1.0e-6, 1.0e-6, 1.0e-15);


    if(print_option>=0)
      err += print_inputs(mat_pvp, sim);
    
    // perform simulations
    double p0 = mat_pvp.p0;
    double h  = poro_visco_plasticity_hardening(p0, &mat_pvp);
    double HardLawJp0Coeff = pow(exp(h), 1.0/3.0);
    
    // deformation gradients
    double Fnp1[DIM_3x3], Fn[DIM_3x3], pFnp1[DIM_3x3], pFn[DIM_3x3], L[DIM_3x3];
    double  F0[DIM_3x3] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double   I[DIM_3x3] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};


    double h_pc_inf  = poro_visco_plasticity_hardening(mat_pvp.pc_inf, &mat_pvp);
    double TMD_0 = exp(h_pc_inf);

    double pJ_0 = TMD_0/sim.TMD;
    double pow_pJ = pow(pJ_0, 1.0/3.0);
    F0[0] = F0[4] = F0[8] =  pow_pJ;    
    
    double pc_n = poro_visco_plasticity_compute_pc(pJ_0,
                                                   (mat_pvp.p0 + mat_pvp.pc_inf)/2.0,
                                                   &mat_pvp,
                                                   &solver_info);

    if(print_option==1){
      cout << "TMD_0 = " << TMD_0   << endl;
      cout << "TMD   = " << sim.TMD << endl;
      cout << "pJ_0  = " << pJ_0    << endl;
      cout << "pc_n  = " << pc_n    << endl;        
    }
    
    double pc_np1 = pc_n;
    
    memcpy(pFn,   F0, sizeof(double)*DIM_3x3);
    memcpy(pFnp1, F0, sizeof(double)*DIM_3x3);
    memcpy( Fn,   F0, sizeof(double)*DIM_3x3);
    memcpy( Fnp1, F0, sizeof(double)*DIM_3x3);    

    char fname[2048];
    sprintf(fname, "%s.txt", sim.file_out); 
    FILE *fp = fopen(fname, "w");

    Tensor<2,3,double *> pF(pFnp1);
    
    struct timeval start, end;
    gettimeofday(&start, NULL);
    
    double L0[3];
    double l[3];
    
    for(int ia=0; ia<3; ia++)
      L0[ia] = Fnp1[ia*3+0] + Fnp1[ia*3+1] + Fnp1[ia*3+2];

    if(print_option==1){  
      cout << "--------------------------------------------" << endl;
      cout << "Simulation results ([time] [pc] [J] [eJ] [pJ] sigma(11~33))" << endl;  
      cout << "--------------------------------------------" << endl;
      cout << "initial geometry (computed from 1x1x1 box) = " << L0[0] << ", " << L0[1] << ", " << L0[2] << endl;
    }
    
    for(int iA=1; iA<=sim.stepno; iA++)
    {
      double t = iA*(sim.dt);
      
      // compute total deformation gradient using velocity gradient
      F_of_t(Fn,Fnp1,L,t,sim);
      
      err += poro_visco_plasticity_integration_algorithm(&mat_pvp, &solver_info, Fnp1, Fn, pFnp1, pFn, &pc_np1, pc_n, sim.dt, sim.implicit);
      err += print_results(fp, Fnp1, pFnp1, pc_np1, t, &mat_pvp, L0, print_option); 

      memcpy(pFn,pFnp1,sizeof(double)*DIM_3x3);
      memcpy( Fn, Fnp1,sizeof(double)*DIM_3x3);      
      pc_n = pc_np1;
      
    }
    fclose(fp);
    
    gettimeofday(&end, NULL);
    double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0 
    + (double)(end.tv_sec - start.tv_sec);
    printf ("Total time: %.4lf s\n", diff);
    
  }
  return err;
}





