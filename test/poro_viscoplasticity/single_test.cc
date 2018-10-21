/// Authors:
/// Sangmin Lee, [1], <slee43@nd.edu>
/// 
/// [1] University of Notre Dame, Notre Dame, IN

#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h> // clock_t, clock, CLOCKS_PER_SEC

#include <ttl/ttl.h>

#include "crystal_plasticity_integration.h"
#include"poro_visco_plasticity.h"

using namespace std;
template<int R=2, int D = 3, class S=double>
using Tensor = ttl::Tensor<R, D, S>;

static constexpr ttl::Index<'i'> i;
static constexpr ttl::Index<'j'> j;
static constexpr ttl::Index<'k'> k;
static constexpr ttl::Index<'l'> l;
    
#define DIM_3x3 9
           
int main(int argc,char *argv[])
{
  double dt = 0.01;
  MaterialPoroViscoPlasticity mat_pvp;
  
  GcmSolverInfo solver_info;
  set_gcm_solver_info(&solver_info, 10, 10, 100, 1.0e-6, 1.0e-6, 1.0e-15);  
  solver_info.debug = true;

  set_properties_poro_visco_plasticity(&mat_pvp, 1.0,   1.1,   0.15,  0.0005, 0.62, 
                                                 0.37, 77.22, 13.01, 15,      0.01, 
                                                 0.2,   5.8,  30,    60,      0.063, 
                                                 0.008, 2,     1,   290);
                                                 
  PvpElasticity elast(&mat_pvp, true);                                                 
                                                 
  double Fnp1[DIM_3x3]  = {9.88552621283890209e-01,0.00000000000000000e+00,0.00000000000000000e+00,
                           0.00000000000000000e+00,9.88552621283890209e-01,0.00000000000000000e+00,
                           0.00000000000000000e+00,0.00000000000000000e+00,9.88552621283890209e-01};
                        
  double Fn[DIM_3x3]    = {9.88584697162311676e-01,0.00000000000000000e+00,0.00000000000000000e+00,
                           0.00000000000000000e+00,9.88584697162311676e-01,0.00000000000000000e+00,
                           0.00000000000000000e+00,0.00000000000000000e+00,9.88584697162311676e-01};
                        
  double pFn[DIM_3x3]   = {9.99066755204096490e-01,0.00000000000000000e+00,0.00000000000000000e+00,
                           0.00000000000000000e+00,9.99066755204096490e-01,0.00000000000000000e+00,
                           0.00000000000000000e+00,0.00000000000000000e+00,9.99066755204096490e-01};
                          
  double pFnp1[DIM_3x3] = {9.99066755204096490e-01,0.00000000000000000e+00,0.00000000000000000e+00,
                           0.00000000000000000e+00,9.99066755204096490e-01,0.00000000000000000e+00,
                           0.00000000000000000e+00,0.00000000000000000e+00,9.99066755204096490e-01};

  double pc_n = 2.66406654285421629;
  double pc_np1 = pc_n;
  
  int err = poro_visco_plasticity_integration_algorithm(&mat_pvp, &solver_info, &elast, Fnp1, Fn, pFnp1, pFn, &pc_np1, pc_n, dt);


  double sol_pFnp1[DIM_3x3] = {9.99047716772131000e-01,0.00000000000000000e+00,0.00000000000000000e+00,
                               0.00000000000000000e+00,9.99047716772131000e-01,0.00000000000000000e+00,
                               0.00000000000000000e+00,0.00000000000000000e+00,9.99047716772131000e-01};

  double sol_pc_np1 = 2.67513433851419169;


  printf("pFnp1(solution)={\n%.17e,%.17e,%.17e,\n%.17e,%.17e,%.17e,\n%.17e,%.17e,%.17e};\n", sol_pFnp1[0], sol_pFnp1[1], sol_pFnp1[2],
                                                                                             sol_pFnp1[3], sol_pFnp1[4], sol_pFnp1[5],
                                                                                             sol_pFnp1[6], sol_pFnp1[7], sol_pFnp1[8]);
  
  printf("pFnp1(computed)={\n%.17e,%.17e,%.17e,\n%.17e,%.17e,%.17e,\n%.17e,%.17e,%.17e};\n", pFnp1[0], pFnp1[1], pFnp1[2],
                                                                                             pFnp1[3], pFnp1[4], pFnp1[5],
                                                                                             pFnp1[6], pFnp1[7], pFnp1[8]);
  printf("pc_np1(sol:computed) = %.17e : %.17e\n", sol_pc_np1, pc_np1);
  return err;
}




