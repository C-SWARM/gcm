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
  PvpElasticity elast;
  Tensor<2,3,double> F, pF;
  double pc;
  PARAMS(){
    elast.construct_elasticity(&mat_pvp,false);
    F = {};
    pF = {};
    pc = 0.0;
  }
  ~PARAMS(){}
};

int print_results(PARAMS &param){
  int err = 0;

  cout << "-----------------------------------------------------------" << endl;
  cout << "Computed stress                                            " << endl;
  cout << "-----------------------------------------------------------" << endl;
  Tensor<2,3,double> eF, pF_I, FI;
  pF_I = ttl::inverse(param.pF);
  eF = param.F(i,k)*pF_I(k,j);
  FI  = ttl::inverse(param.F);

  double eJ = ttl::det(eF);
  double pJ = ttl::det(param.pF);
  double  J = ttl::det(param.F);

  Tensor<2,3,double> S, eS, sigma;
  param.elast.set_internal_variable(param.pc);
  param.elast.update_elasticity(eS.data, NULL, eF.data, false);
  sigma(i,l) = eF(i,j)*eS(j,k)*eF(l,k)/eJ;
  
  S(i,j) = J*FI(i,k)*sigma(k,l)*FI(j,l);
  
  cout<<"pc = " << param.pc << ", J = "<< J << ", eJ = " << eJ << ", pJ = " << pJ << endl;
  cout<< "sigma = {";

  for(int ia=0; ia<DIM_3x3-1; ++ia)
    cout << sigma.get(ia) << ", ";
  
  cout << sigma.get(DIM_3x3-1) << "};" << endl;
  
  cout<< "PKII = {";

  for(int ia=0; ia<DIM_3x3-1; ++ia)
    cout << S.get(ia) << ", ";
  
  cout << S.get(DIM_3x3-1) << "};" << endl;  
      
  return err;
}
                  
/// read input file
///
/// Input file is formated as:
/// # is comments
/// # is comments
/// #-------------------------------------------------------------------------
/// #  M   alpha m_d  m_v  gamma0 a1   a2   Lambda1 Lambda2  c_inf Gamma B   pcb d_m mu_0 mu_1 p0    kappa n g0 pc_inf
/// #-------------------------------------------------------------------------
///    1.0 1.1   0.15 0.15 0.0005 0.62 0.37 77.22   13.01    15    0.01  0.2 5.8 1   30   60   0.063 0.008 2 1  290
/// # pc
/// 2.67513433851419169       
///
/// # total deformation gradient
/// 9.88552621283890209e-01 0.00000000000000000e+00 0.00000000000000000e+00
/// 0.00000000000000000e+00 9.88552621283890209e-01 0.00000000000000000e+00
/// 0.00000000000000000e+00 0.00000000000000000e+00 9.88552621283890209e-01
///
/// # plastic deformation gradient
/// 9.99047716772131000e-01 0.00000000000000000e+00 0.00000000000000000e+00
/// 0.00000000000000000e+00 9.99047716772131000e-01 0.00000000000000000e+00
/// 0.00000000000000000e+00 0.00000000000000000e+00 9.99047716772131000e-01
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
  fscanf(fp, "%lf", &param.pc);

  err += goto_valid_line(fp);
  for(int ia=0; ia<DIM_3x3; ++ia)
    fscanf(fp, "%lf", param.F.data + ia);
    
  err += goto_valid_line(fp);
  for(int ia=0; ia<DIM_3x3; ++ia)
    fscanf(fp, "%lf", param.pF.data + ia);
    
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
  
  cout << "F = {";
  for(int ia=0; ia<DIM_3x3-1; ++ia)
    cout << param.F.data[ia] << ", ";

  cout << param.F.data[DIM_3x3-1] << "};" << endl;
    
  cout << "pF = {";
  for(int ia=0; ia<DIM_3x3-1; ++ia)
    cout << param.pF.data[ia] << ", ";

  cout << param.pF.data[DIM_3x3-1] << "};" << endl;

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





