/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
/// 
/// [1] University of Notre Dame, Notre Dame, IN

#include <iostream>
#include "input_file_read_handles.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "continuum_damage_model.h"
#include "crystal_plasticity_integration.h"
#include <sys/time.h>
#include <string.h>
#include <ttl/ttl.h>

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

/// print damage model results
///
/// \param[in]  fp_S          file pointer for writing simulation results(Stress)
/// \param[in]  fp_F          file pointer for writing simulation results(Deformation gradient)
/// \param[in]  F             total deformation gradient
/// \param[in]  w             damage
/// \param[in]  H             dGdY from the consistency condition
/// \param[in]  is_it_damaged if 1 : damage was evolved
///                             0 : damage wasn't evolved
/// \param[in]  dt            time step size
/// \param[in]  t             time
/// \param[in]  mat_d         material property object for damage model
/// \param[out] elast         elasticity object. elast.S will be updated
/// \param[in]  print_option  if 1: display computed values
void print_results(FILE *fp_S,
                   FILE *fp_F,
                   ttl::Tensor<2,3,double> &F,
                   const double dw,
                   const double dH,
                   const double vw,
                   const double vH,
                   const int is_it_damaged_d,
                   const int is_it_damaged_v,
                   const double dt,
                   const double t,
                   MATERIAL_CONTINUUM_DAMAGE &mat_d,
                   HyperElasticity &elast,
                   const int print_option){

  bool construct_ET = true;
  update_split_damage_elasticity(&mat_d,&elast,dw,vw,dH,vH,
                                 is_it_damaged_d,is_it_damaged_v,
                                 dt, F.data, construct_ET);
                                        

  ttl::Tensor<2,3,double> sigma;
  elast.compute_Cauchy(sigma.data,F.data);
    
  double  J = ttl::det(F);

  if(print_option==1)        
    printf("%.17e %.17e %.17e %.17e %.17e %.17e\n", t, dw, vw, elast.S[0],elast.S[4],elast.S[8]);

  fprintf(fp_S, "%.17e %.17e %.17e %.17e", t, dw, vw, J);

  for(int ia=0; ia<DIM_3x3; ia++)
    fprintf(fp_S, "%.17e ", sigma.get(ia));

  fprintf(fp_S, "\n");
      
  for(int ia=0; ia<DIM_3x3; ia++)
    fprintf(fp_F, "%.17e ", F.get(ia));
    
  fprintf(fp_F, "\n");
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
/// #  E     nu    devPotFlag volPotFlag mu  ome_max[0,1) p1  p2  Yin
/// #  alpha_dev beta_dev alpha_vol beta_vol
/// #-------------------------------------------------------------------------
///    800   0.34  1          2          100 1.0          8.0 2.5 0.15
///    1 1 1 1
/// #-------------------------------------------------------------------------
/// # Analysis name
/// #-------------------------------------------------------------------------
///   compression
/// #
/// # time steps
/// # dt   number_of_time_steps
/// #-------------------------------------------------------------------------
///   1.0  100
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
/// \param[in, out] &mat_e   material property object for elasticity
/// \param[in, out] &mat_d   material property object for damage model
/// \param[in, out] &sim     simulation parameter object such as numer of time step, time step size
/// \return         non-zero on interal error
int read_input_file(const char *input_file,
                    MATERIAL_ELASTICITY &mat_e,
                    MATERIAL_CONTINUUM_DAMAGE &mat_d,
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

  double E, nu, mu, ome_max, p1, p2, Yin;
  double alpha_dev, beta_dev, alpha_vol, beta_vol;
  int devPotFlag, volPotFlag;
  
  err += goto_valid_line(fp_sim);
      
  // read material parameters 
  fscanf(fp_sim, "%lf %lf", &E, &nu);
  fscanf(fp_sim, "%d %d", &devPotFlag, &volPotFlag);  
  fscanf(fp_sim, "%lf %lf %lf %lf %lf", &mu, &ome_max, &p1, &p2, &Yin);
  fscanf(fp_sim, "%lf %lf %lf %lf", &alpha_dev, &beta_dev, &alpha_vol, &beta_vol);
    
  set_properties_using_E_and_nu(&mat_e, E, nu);
  mat_e.devPotFlag = devPotFlag;
  mat_e.volPotFlag = volPotFlag;
    
  set_split_damage_parameters(&mat_d, p1, p2, Yin, mu, ome_max,
                              alpha_dev, beta_dev, alpha_vol, beta_vol);  
  
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
/// \param[in] &mat_e   material property object for elasticity
/// \param[in] &mat_d   material property object for damage model
/// \param[in] &sim     simulation parameter object such as numer of time step, time step size
/// \return         non-zero on interal error 
int print_inputs(MATERIAL_ELASTICITY &mat_e,
                 MATERIAL_CONTINUUM_DAMAGE &mat_d,
                 SIM_PARAMS &sim)
{
  int err = 0;

  cout << "--------------------------------------------" << endl;
  cout << "Numerical parameters" << endl;  
  cout << "--------------------------------------------" << endl;            
  cout << "output file name\t: "     << sim.file_out  << endl;
  cout << "dt\t\t\t: "               << sim.dt        << endl;
  cout << "number of time steps\t: " << sim.stepno    << endl;
  cout << "loading velocity \t: "    << sim.velocity  << endl;
  cout << "loading dimension \t: "   << sim.dim       << endl;  
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
  print_material_property_elasticity(&mat_e);
  print_material_property_damage_model(&mat_d);
  

  return err;
} 

int main(int argc,char *argv[])
{
  int err = 0;
    
  char fn_sim[1024];
  int print_option = 0;

  if(argc<2)
  {
    cout << "Usage:./test_damage_model [FILE_SIM] [print_option]\n";
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
      
      
    MATERIAL_ELASTICITY mat_e;
    MATERIAL_CONTINUUM_DAMAGE mat_d;
    HyperElasticity elast;
    SIM_PARAMS sim;

    err += read_input_file(fn_sim, mat_e, mat_d, sim);
    elast.construct_elasticity(&mat_e, true);
    
    if(print_option>=0)
      err += print_inputs(mat_e, mat_d, sim);

    // deformation gradients
    ttl::Tensor<2,3,double> L  = {};
    ttl::Tensor<2,3,double> F  = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
    ttl::Tensor<2,3,double> Fn = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};

    char fname_S[2048], fname_F[2048];
    sprintf(fname_S, "%s.sigma.txt", sim.file_out); 
    sprintf(fname_F, "%s.F.txt", sim.file_out); 
        
    FILE *fp_S = fopen(fname_S, "w");
    FILE *fp_F = fopen(fname_F, "w");    

    if(print_option==1)
    {  
      cout << "--------------------------------------------" << endl;
      cout << "Simulation results ([time] [wd] [wv] sigma(11~33))" << endl;  
      cout << "--------------------------------------------" << endl;      
    }
        
    struct timeval start, end;
    gettimeofday(&start, NULL);
  
    double dw, dX, dH, dwn, dXn;
    dw = dwn = dX = dXn = dH = 0.0;
    double vw, vX, vH, vwn, vXn;
    vw = vwn = vX = vXn = vH = 0.0;    

    int is_it_damaged_d = 0;
    int is_it_damaged_v = 0;    

    for(int iA=1; iA<=sim.stepno; iA++){
      double t = iA*(sim.dt);
      // compute total deformation gradient using velocity gradient
      F_of_t(Fn.data,F.data,L.data,t,sim);      
      
      err += continuum_split_damage_integration_alg(&mat_d,&elast,
                                                    &dw,&vw,&dX,&vX,&dH,&vH,
                                                    &is_it_damaged_d,&is_it_damaged_v,                                     
                                                    dwn,vwn,dXn,vXn,sim.dt,F.data); 
                                           
      dwn = dw;
      dXn = dX;
      vwn = vw;
      vXn = vX;      
      Fn(i,j) = F(i,j);
      is_it_damaged_d = is_it_damaged_v = 0;

      print_results(fp_S, fp_F, F, dw, dH, vw, vH, is_it_damaged_d, is_it_damaged_v,
                    sim.dt, t, mat_d, elast, print_option);
    }
    fclose(fp_S);
    fclose(fp_F);    
        
    gettimeofday(&end, NULL);
    double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0 
    + (double)(end.tv_sec - start.tv_sec);
    printf ("Total time: %.4lf s\n", diff);
    
  }
}
