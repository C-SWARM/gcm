#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "crystal_plasticity_integration.h"
#include "flowlaw.h"
#include <sys/time.h>

void test_crystal_plasticity_single_crystal(void)
{
  //double lame1 = 75600.0;
  //double lame2     = 26100.0;
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
  double tol_hardening  = 1.0e-4;
  double tol_M          = 1.0e-4;
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
  CRYSTAL_PLASTICITY_SOLVER_INFO solver_info;
  set_crystal_plasticity_solver_info(&solver_info,max_itr_stag,
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
  Tensor<2> pFn,pFnp1,pFnp1_I,eFnp1,Fn,Fnp1,L={},sigma,PK2dev,sigma_dev;
  
  pFn   = ttl::identity(i,j);
  pFnp1 = ttl::identity(i,j);
  eFnp1 = ttl::identity(i,j);
  Fn    = ttl::identity(i,j);
  Fnp1  = ttl::identity(i,j);
      
  double g_n,g_np1;
  g_n = g_np1 = mat_p.g0;
  
  double dt = 0.001;
  
  double d = 1.0;
  // set velocity gradient  
  L[0][0] = -d;
  L[1][1] = L[2][2] = d/2;  
  
  // start integration  
  TensorA<2> PK2(elast.S);
  
  FILE *fp = fopen("single_crystal_results.txt", "w");
  
  for(int a = 1; a<=1000; a++)
  {
    double lambda = 0.0;
    double t = a*dt;
    
    // compute total deformation gradient using velocity gradient
    Fnp1_Implicit(Fnp1.data, Fn.data, L.data, dt); 
    
    staggered_Newton_Rapson(pFnp1.data, &g_np1, &lambda, 
                            pFn.data, Fn.data,Fnp1.data, 
                            g_n, dt, &mat, &elast, &solver_info);
                            
    pFn = pFnp1(i,j);
    Fn  = Fnp1(i,j);
    inv(pFnp1,pFnp1_I);
    eFnp1 = Fnp1(i,k)*pFnp1_I(k,j);
       
    g_n = g_np1;
    
    
    // print result at time t
    double sigma_eff;
    double PK2_eff;
    double det_pF = ttl::det(pFnp1);
        
    elast.update_elasticity(&elast,eFnp1.data,0);
    elast.compute_PK2_eff(&elast,&PK2_eff);
    elast.compute_Cauchy_eff(&elast,&sigma_eff,eFnp1.data);   

    fprintf(fp, "%e %e %e %e %e %e\n",t,sigma_eff,PK2_eff, g_np1, det_pF, PK2[0][0]);
  }    
  
  fclose(fp);  
  destruct_elasticity(&elast);
  destruct_slip_system(&slip);
}


int main(int argc,char *argv[])
{
  int err = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  test_crystal_plasticity_single_crystal();
  
  gettimeofday(&end, NULL);
  double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0 
  + (double)(end.tv_sec - start.tv_sec);
  printf ("Total time: %.4lf s\n", diff);
  
  return err;
}
