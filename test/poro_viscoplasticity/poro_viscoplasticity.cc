//
//  main.cpp
//  ttl-learning
//
//  Created by alberto salvadori on 12/8/16.
//  Copyright © 2016 alberto salvadori. All rights reserved.
//

// include

#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h> // clock_t, clock, CLOCKS_PER_SEC

#include <ttl/ttl.h>

#include "KMS-IJSS2017.h"
#include "ttl-tools.h"
#include "input_file_read_handles.h"
#include "hyperelasticity.h"
#include "crystal_plasticity_integration.h"
#include "pvp_interface.h"

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
  PARAM_pc_0,       // initial pc
  PARAM_pJ,         // initial plastic deformation
  PARAM_NO
};

int compute_stess(double *sigma_in,
                  double *F_in,
                  double *pF_in,
                  double *S_in,
                  double pc,
                  KMS_IJSS2017_Parameters &mat_pvp)
{
  int err = 0;
  Tensor<2,3,double*> S(S_in), F(F_in), pF(pF_in),sigma(sigma_in);
  Tensor<2,3,double> eF, pF_I;
  pF_I = ttl::inverse(pF);
  eF = F(i,k)*pF_I(k,j);
  double eJ = ttl::det(eF);
  
  err += pvp_intf_update_elasticity(eF.data,pc,S_in,NULL,&mat_pvp,0);

  sigma(i,l) = eF(i,j)*S(j,k)*eF(l,k)/eJ;                                        
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
/// #-------------------------------------------------------------------------
/// #  yf_M  yf_alpha  flr_m flr_gamma0  hr_a1 hr_a2 hr_Lambda1 hr_Lambda2  c_inf c_Gamma d_B d_pcb mu_0 mu_1 K_p0  K_kappa pl_n cf_g0 cf_pcinf
/// #-------------------------------------------------------------------------
///    1.0   1.1       0.15  0.0005      0.62  0.37  77.22      13.01       15    0.01    0.2 5.8   30   60   0.063 0.008   2    1     290
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
/// # 0: F_of_t(L1,L2,L3, ...) : L = n number of velocity gradient with t(end) to be transient when t reach t_end
/// # 1: F_of_t(v,dim)         : v = velocity (displacement), dim = 1: uniaxial displacement
/// #                                                               2: biaxial displacement
/// #                                                               3: triaxial displacement
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
                    KMS_IJSS2017_Parameters &mat_pvp,
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
  
  // read material parameters 
  double param[PARAM_NO];
  for(int ia=0; ia<PARAM_NO; ia++)
    fscanf(fp_sim, "%lf", param+ia);
  
  mat_pvp.set_parameters(param[PARAM_yf_M],
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
                         param[PARAM_cf_pcinf],
                         true);
  
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
int print_inputs(KMS_IJSS2017_Parameters &mat_pvp,
                 SIM_PARAMS &sim)
{
  int err = 0;

  cout << "--------------------------------------------" << endl;
  cout << "Numerical parameters" << endl;  
  cout << "--------------------------------------------" << endl;
  cout << "output file name\t: "      << sim.file_out  << endl;
  cout << "dt\t\t\t: "               << sim.dt        << endl;
  cout << "number of time steps\t: " << sim.stepno    << endl;
  cout << "integration type\t: "     << sim.intg_type << endl;
  cout << "loading velocity \t: "  << sim.velocity  << endl;
  cout << "loading dimension \t: " << sim.dim       << endl;  
  cout << "number of Ls\t\t: "       << sim.Lno       << endl;
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
  std::string s;
  mat_pvp.AsAString(s);
    
  cout << s;

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
    
    KMS_IJSS2017_Parameters mat_pvp;
    SIM_PARAMS sim;

    err += read_input_file(fn_sim, mat_pvp, sim);

    if(print_option>=0)
      err += print_inputs(mat_pvp, sim);
    
    // perform simulations
    double p0 = mat_pvp.get_K_p0();
    double h  = pvp_intf_hardening_law(p0, &mat_pvp);
    double HardLawJp0Coeff = pow(exp(h), 1.0/3.0);
    
    // deformation gradients
    //double L[DIM_3x3];
    // stress
    double PKII[DIM_3x3], sigma[DIM_3x3]; 
    double  F0[DIM_3x3] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    //double   I[DIM_3x3] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};                  

    F0[0] = F0[4] = F0[8] =  HardLawJp0Coeff;
                 
    double pc_n, pc_np1;
    pc_n = pc_np1 = p0;



double dt = 1.666667e-10;
pc_n = 6.454024e+0;
double Fnp1[DIM_3x3] = {1.01683031257434564e+00,1.22170066969445532e-03,7.16606936255713856e-03,9.98743992990401563e-04,1.00897793558907489e+00,8.80184670223167331e-04,-2.12552201971938359e-03,8.66326795157812413e-03,6.41702793772817914e-01};
double Fn[DIM_3x3] = {1.01683031213029884e+00,1.22170069481151999e-03,7.16606977226791862e-03,9.98743948766325773e-04,1.00897793587215689e+00,8.80184645444615941e-04,-2.15893815462742111e-03,8.70965779481824853e-03,6.41729621578168485e-01};
double pFnp1[DIM_3x3] = {2.17018427086157129e+00,7.80158597775989365e-02,8.73913759600668527e-02,7.99005570420286554e-02,1.39707416846036181e+00,6.58599687252051091e-02,1.19813588812333047e-01,9.76534909234337001e-02,3.75631034347642390e-01};
double pFn[DIM_3x3] = {9.99358697027337595e-01,-4.06813610025157893e-04,3.13044484217239218e-03,1.06843278872218849e-03,1.00150237906777306e+00,7.37428680410462281e-04,-3.17738984603657985e-03,9.19446648932062470e-05,6.62367297566719193e-01};



//double dt = 5.709247e-15;
//pc_n = 6.454024e+01;
//double Fnp1[DIM_3x3] = {1.01679482240293151e+00,1.22703867437358102e-03,7.24793264027288626e-03,1.03924000226460200e-03,1.00902010024838673e+00,9.18399776197452097e-04,-2.17844258198152119e-03,8.71286857042430363e-03,6.41722886484102806e-01};
//double Fn[DIM_3x3] = {1.01679482240293151e+00,1.22703867437358102e-03,7.24793264027288626e-03,1.03924000226460200e-03,1.00902010024838673e+00,9.18399776197452097e-04,-2.17844372710107592e-03,8.71287016012045612e-03,6.41722887403513020e-01};
//double pFnp1[DIM_3x3] = {2.17014948370544181e+00,7.80344115859264165e-02,8.75707762513473764e-02,7.99550616215682075e-02,1.39714247676462522e+00,6.59197761424780010e-02,1.19803929731099129e-01,9.76631848603330965e-02,3.75639024227740981e-01};
//double pFn[DIM_3x3] = {9.99323873994986900e-01,-4.01651947110411576e-04,3.21088233952578853e-03,1.10863035910188038e-03,1.00154423083781485e+00,7.75368344947294924e-04,-3.19784741297175030e-03,9.48836812957112863e-05,6.62359934012867413e-01};    
    char fname[1024];
    sprintf(fname, "%s.txt", sim.file_out); 
    FILE *fp = fopen(fname, "w");
    
    if(print_option==1)
    {  
      cout << "--------------------------------------------" << endl;
      cout << "Simulation results ([time] [pC] [pJ] sigma(11))" << endl;  
      cout << "--------------------------------------------" << endl;      
    }
    
    Tensor<2,3,double *> pF(pFnp1);




    
    for(int iA=1; iA<=1; iA++)
    {
      double t = iA*(sim.dt);
  
      err += pvp_intf_perform_integration_algorithm(Fnp1,Fn,pFnp1,pFn,&pc_np1,pc_n,&mat_pvp,dt);
      double pJ = ttl::det(pF);

      err += compute_stess(sigma,Fnp1,pFnp1,PKII,pc_np1,mat_pvp);
                  

      if(print_option==1)        
        printf("%.17e, %.17e %.17e %.17e\n", dt, pc_np1, pJ, sigma[0]);
        
      //print_T2(pFnp1, "pFnp1");
      //print_T2(pFn,   "pFn");  

      fprintf(fp, "%.17e, %.17e %.17e %.17e\n", t, pc_np1, pJ, sigma[0]); 
      memcpy(pFn,pFnp1,sizeof(double)*DIM_3x3);
      memcpy( Fn, Fnp1,sizeof(double)*DIM_3x3);
      pc_n = pc_np1;
    }
    fclose(fp);
  }
  return err;
}





