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

using namespace std;
template<int R=2, int D = 3, class S=double>
using Tensor = ttl::Tensor<R, D, S>;

static constexpr ttl::Index<'i'> i;
static constexpr ttl::Index<'j'> j;
static constexpr ttl::Index<'k'> k;
static constexpr ttl::Index<'l'> l;
    
#define DIM_3x3   9
#define PARAM_NO 19

class PARAMS
{
public:
  MaterialPoroViscoPlasticity mat_pvp;
  double pJ;
  PARAMS(){
    pJ = 1.0;
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
    
  double pc = poro_visco_plasticity_compute_pc(param.pJ,
                                               param.mat_pvp.p0,
                                               &param.mat_pvp,
                                               &solver_info);

  cout<<"pc = " << pc << ", pJ = "<< param.pJ << endl;
  return err;
}
                  
/// read input file
///
/// Input file is formated as:
/// # is comments
/// #-------------------------------------------------------------------------
/// #  yf_M  yf_alpha  flr_m flr_gamma0  hr_a1 hr_a2 hr_Lambda1 hr_Lambda2  c_inf c_Gamma d_B d_pcb mu_0 mu_1 K_p0  K_kappa pl_n cf_g0 cf_pcinf
/// #-------------------------------------------------------------------------
///    1.0   1.1       0.15  0.0005      0.62  0.37  77.22      13.01       15    0.01    0.2 5.8   30   60   0.063 0.008   2    1     290
/// # pJ
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
  double pvp_param[PARAM_NO];
  for(int ia=0; ia<PARAM_NO; ia++)
    fscanf(fp, "%lf", pvp_param+ia);
    
  set_properties_poro_visco_plasticity(&param.mat_pvp,
                                       pvp_param[ 0], pvp_param[ 1], pvp_param[ 2], pvp_param[ 3],
                                       pvp_param[ 4], pvp_param[ 5], pvp_param[ 6], pvp_param[ 7],
                                       pvp_param[ 8], pvp_param[ 9], pvp_param[10], pvp_param[11],
                                       pvp_param[12], pvp_param[13], pvp_param[14], pvp_param[15],
                                       pvp_param[16], pvp_param[17], pvp_param[18]);
  
  err += goto_valid_line(fp);
  fscanf(fp, "%lf", &param.pJ);
    
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
  
  cout << "pJ = " << param.pJ << endl;
  return err;
}                 

int main(int argc,char *argv[])
{
  int err = 0;
    
  char fn[1024];

  if(argc<2)
  {
    cout << "Usage: ./compute_stress [FILE_PARAM]\n";
    cout << "\t[FILE_PARAM]\t: file path for parameters\n";
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





