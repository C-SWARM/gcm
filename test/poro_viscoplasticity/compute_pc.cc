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
#include "poro_visco_plasticity.h"
#include "read_params.h" // local head providing reading parameters

using namespace std;
template<int R=2, int D = 3, class S=double>
using Tensor = ttl::Tensor<R, D, S>;

static constexpr ttl::Index<'i'> i;
static constexpr ttl::Index<'j'> j;
static constexpr ttl::Index<'k'> k;
static constexpr ttl::Index<'l'> l;
    
#define DIM_3x3   9

class PARAMS
{
public:
  MaterialPoroViscoPlasticity mat_pvp;
  double TMD;
  PARAMS(){
    TMD = 1.0;
  }
  ~PARAMS(){}
};

int print_results(PARAMS &param){
  int err = 0;

  cout << "-----------------------------------------------------------" << endl;
  cout << "Computed pc                                                " << endl;
  cout << "-----------------------------------------------------------" << endl;

  GcmSolverInfo solver_info;
  set_gcm_solver_info(&solver_info, 10, 10, 100, 1.0e-6, 1.0e-6, 1.0e-15);
  solver_info.debug = true;


  double h_pc_inf  = poro_visco_plasticity_hardening(param.mat_pvp.pc_inf, &param.mat_pvp);
  double TMD_0 = exp(h_pc_inf);

  double pJ = TMD_0/param.TMD;
  double pc = poro_visco_plasticity_compute_pc(pJ,
                                              (param.mat_pvp.p0 + param.mat_pvp.pc_inf)/2.0,
                                              &param.mat_pvp,
                                              &solver_info);
    
  cout << "TMD = "<< param.TMD << ", TMD0 = " << TMD_0 << ", pJ = "<< pJ <<": pc = " << pc << endl;
  return err;
}
                  
/// read input file
///
/// Input file is formated as:
/// # is comments
/// #-------------------------------------------------------------------------
/// #  M   alpha m_d  m_v  gamma0 a1   a2   Lambda1 Lambda2  c_inf Gamma B   pcb d_m mu_0 mu_1 p0    kappa n g0 pc_inf
/// #-------------------------------------------------------------------------
///    1.0 1.1   0.15 0.15 0.0005 0.62 0.37 77.22   13.01    15    0.01  0.2 5.8 1   30   60   0.063 0.008 2 1  290
/// # TMD
/// 0.9
/// \return non-zero on interal error
int read_input_file(const char *input_file,
                    PARAMS &param)
{
  int err = 0;
  
  FILE *fp = NULL;
  fp = fopen(input_file, "r");

  // read material and simulation parameters
  // if fail to read the file, exit    
  if(fp==NULL)
  { 
    cout << "fail to read ["<< input_file << "]\n";
    exit(-1);      
  }  

  err += goto_valid_line(fp);

  // read material parameters
  read_pvp_material_properties(fp, param.mat_pvp);
  
  err += goto_valid_line(fp);
  fscanf(fp, "%lf", &param.TMD);
    
  fclose(fp);
  
  return err;
}                    

/// print numerical and material parameters
/// 
/// \param[in, out] &mat_pvp material property object for porovisco-plasticity
/// \param[in, out] &param   imulation parameter object such as numer of time step, time step size
/// \return         non-zero on interal error 
int print_inputs(PARAMS &param)
{
  int err = 0;

  // print Material parameters  
  print_material_property_poro_visco_plasticity(&param.mat_pvp);
  
  cout << "TMD = " << param.TMD << endl;
  return err;
}                 

int main(int argc,char *argv[])
{
  int err = 0;
    
  char fn[1024];

  if(argc<2)
  {
    cout << "Usage: ./compute_stress [FILE_PARAM]" << endl;
    cout << "\t[FILE_PARAM]\t: file path for parameters" << endl;
    cout << "\t\t\t#Input file is formated as:" << endl;
    cout << "\t\t\t# is comments" << endl;
    cout << "\t\t\t#-------------------------------------------------------------------------" << endl;
    cout << "\t\t\t#  M    alpha  m_d  m_v  gamma0 a1    a2   Lambda1 Lambda2  c_inf Gamma  B   pcb d_m mu_0 mu_1 K_p0  K_kappa n g0 pc_inf" << endl;
    cout << "\t\t\t#-------------------------------------------------------------------------" << endl;
    cout << "\t\t\t   1.0  1.1    0.15 0.15 0.0005 0.62  0.37 77.22   13.01    15    0.01   0.2 5.8  1  30   60   0.063 0.008   2 1  290" << endl;
    cout << "\t\t\t# TMD" << endl;
    cout << "\t\t\t0.9" << endl;

    exit(-1);
  }
  else
  {    
    // read commend line arguments
    sscanf(argv[1], "%s", fn);
        
    PARAMS param;
    err += read_input_file(fn, param);
    err += print_inputs(param);
    err += print_results(param);
  }
  return err;
}





