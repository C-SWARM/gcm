#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "crystal_plasticity_integration.h"
#include "flowlaw.h"

template<class T1, class T2>
void test_crystal_plasticity_single_crystal(T1 &Fnp1, T2 &hFnp1, const char *filename)
{
  double E = 70.0e+3;
  double nu = 0.25;
  
  double gamma_dot_0 = 1.0;
  double gamma_dot_s = 50.0e+9;
  double m           = 0.05;  
  double g0          = 210.0;
  double G0          = 200.0;
  double gs_0        = 330.0;
  double w           = 0.005;
  
  int max_itr_stag      = 100;
  int max_itr_hardening = 5;
  int max_itr_M         = 100;
  double tol_hardening  = 1.0e-6;
  double tol_M          = 1.0e-6;
  double computer_zero  = 1.0e-15;

  // create material properties: Elasticity
  MATERIAL_ELASTICITY mat_e;  
  set_properties_using_E_and_nu(&mat_e,E,nu); 
  // or you can use : set_properties_using_Lame_constants(&mat_e,lame1,lame2);
  //print_material_property_elasticity(&mat_e); // <= this is optional
  
  // create slip system : 0 for FCC
  // it is needed to setup plasticity
  SLIP_SYSTEM slip;
  construct_slip_system(&slip,0); 
  
  // create material properties: Plasticity
  MATERIAL_CRYSTAL_PLASTICITY mat_p;
  set_properties_crystal_plasticity(&mat_p,&slip,gamma_dot_0,gamma_dot_s, 
                                     m,g0,G0,gs_0,w);
  //print_material_property_crystal_plasticity(&mat_p);  // <= this is optional 

  // create material plasticity: it needs material properties for elasticity and plasticity
  MATERIAL_CONSTITUTIVE_MODEL mat;
  set_properties_constitutive_model(&mat,&mat_e,&mat_p);
  
  // create solver info: criteria for numerical iterations
  GcmSolverInfo solver_info;
  set_gcm_solver_info(&solver_info,max_itr_stag,
                      max_itr_hardening,
                      max_itr_M,
                      tol_hardening,
                      tol_M,
                      computer_zero);  
  //print_crystal_plasticity_solver_info(&solver_info); // <= this is optional
  
  // create elasticity object for integration
  // this creates memory for stress and elasticity tensor s.t. requres destructor
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);  

  // set variables for integration
  Tensor<2> pFnp1,xFn;
  pFnp1 = ttl::identity(i,j);
  xFn   = ttl::identity(i,j);
      
  double g_np1;
  double dt     = 0.1;
  double lambda = 0.0;
  
  GcmCpIntegrator integrator;
  
  integrator.mat         = &mat;
  integrator.elasticity  = &elast;
  integrator.solver_info = &solver_info;  
  
  integrator.set_tensors(Fnp1.data,
                         xFn.data,
                         xFn.data,
                         pFnp1.data,
                         xFn.data,
                         hFnp1.data,
                         xFn.data);

  integrator.gnp1   = &g_np1;
  integrator.lambda = &lambda;
  integrator.gn     = g_np1 = mat_p.g0;
  
  // perform integration algorithm for the crystal plasticity
  integrator.run_integration_algorithm(dt, dt);

  printf("pFnp1 = [\n%e %e %e \n%e %e %e\n%e %e %e]\n",
          pFnp1[0][0], pFnp1[0][1], pFnp1[0][2],
          pFnp1[1][0], pFnp1[1][1], pFnp1[1][2],
          pFnp1[2][0], pFnp1[2][1], pFnp1[2][2]);
              
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "%e %e %e\n%e %e %e\n%e %e %e", pFnp1[0][0], pFnp1[0][1], pFnp1[0][2],
                                              pFnp1[1][0], pFnp1[1][1], pFnp1[1][2],
                                              pFnp1[2][0], pFnp1[2][1], pFnp1[2][2]);
 
  fclose(fp);  
  
  destruct_elasticity(&elast);
  destruct_slip_system(&slip);
}


int main(int argc,char *argv[])
{
  int err = 0;
  
  Tensor<2> Fnp1, hFnp1;
  Fnp1  = ttl::identity(i,j);
  hFnp1 = ttl::identity(i,j);
  
  Fnp1[0][1] = Fnp1[0][2] = 0.01;
  test_crystal_plasticity_single_crystal(Fnp1, hFnp1, "without_thermal.txt");

  hFnp1[0][0] = hFnp1[1][1] = hFnp1[2][2] = 2.0;
   Fnp1[0][0] =  Fnp1[1][1] =  Fnp1[2][2] = 2.0;
   Fnp1[0][1] =  Fnp1[0][2] =  0.02;

  test_crystal_plasticity_single_crystal(Fnp1, hFnp1, "with_thermal.txt");
  
  return err;
}
