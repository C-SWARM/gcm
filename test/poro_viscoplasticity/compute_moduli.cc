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
#include "read_params.h" // local head providing reading parameters

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
  char   file_out[1024]; /// output file name

  SIM_PARAMS(){}
  ~SIM_PARAMS(){}
};
              
/// read input file
///
/// Input file is formated as:
/// # is comments
/// 1 # if 1: implicit 
/// #      0: explicit
/// # is comments
/// #-------------------------------------------------------------------------
/// #  M   alpha m_d  m_v  gamma0 a1   a2   Lambda1 Lambda2  c_inf Gamma B   pcb d_m mu_0 mu_1 p0    kappa n g0 pc_inf
/// #-------------------------------------------------------------------------
///    1.0 1.1   0.15 0.15 0.0005 0.62 0.37 77.22   13.01    15    0.01  0.2 5.8 1   30   60   0.063 0.008 2 1  290
/// # TMD
/// 0.9
/// #-------------------------------------------------------------------------
/// # Analysis name
/// #-------------------------------------------------------------------------
///   iso_compaction
///
///
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
      
  err += goto_valid_line(fp_sim);
  
  // read material parameters
  read_pvp_material_properties(fp_sim, mat_pvp);  

  err += goto_valid_line(fp_sim);
  fscanf(fp_sim, "%s", sim.file_out);

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
  cout << "Numerical parameters"                         << endl;  
  cout << "--------------------------------------------" << endl;            
  cout << "output file name\t: "     << sim.file_out     << endl;  
  cout << "--------------------------------------------" << endl;
  cout << "Material parameters"                          << endl;  
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

    double h_pc_inf  = poro_visco_plasticity_hardening(mat_pvp.pc_inf, &mat_pvp);
    double TMD_0 = exp(h_pc_inf);
    
    if(print_option>=0)
      cout << "TMD_0 = " << TMD_0   << endl;

    double pc_0   = mat_pvp.p0;
    double pc_inf = mat_pvp.pc_inf;
    
    char fname[2048];
    sprintf(fname, "%s.moduli", sim.file_out); 
    FILE *fp = fopen(fname, "w");

    if(print_option>=0){  
      cout << "--------------------------------------------" << endl;
      cout << "Computing moduli: ([TMD] [E] [G] [kappa] [nu])" << endl;  
      cout << "--------------------------------------------" << endl;      
    }
            
    struct timeval start, end;
    gettimeofday(&start, NULL);
    
    double d_pc = (pc_inf - pc_0)/99.0;
    PvpElasticity elast(&mat_pvp, false);
    
    for(int iA=0; iA<100; iA++)
    {
      double pc = pc_0 + d_pc*iA;
        
      double h_pc  = poro_visco_plasticity_hardening(pc, &mat_pvp);
      double pJ = exp(h_pc);
      double TMD = TMD_0/pJ;
      
      elast.set_internal_variable(pc);
      double dudj = 0.0;
      elast.compute_dudj(&dudj, 0.999);
      
      double kappa = -dudj/0.001;
      double mu = poro_visco_plasticity_compute_shear_modulus(&mat_pvp, pc);
      double E = 9.0*kappa*mu/(3.0*kappa+mu);
      
      double nu    = E/2.0/mu - 1.0;
            
      fprintf(fp, "%.17e %.17e %.17e %.17e %.17e\n", TMD, E, mu, kappa, nu);
      if(print_option==1)
        printf("%.17e %.17e %.17e %.17e %.17e\n", TMD, E, mu, kappa, nu);
    }
    fclose(fp);
    
    if(print_option>=0){
      gettimeofday(&end, NULL);
      double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0
                  + (double)(end.tv_sec - start.tv_sec);
      printf ("Total time: %.4lf s\n", diff);
    }    
  }
  return err;
}





