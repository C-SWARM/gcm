#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "crystal_plasticity_integration.h"
#include "flowlaw.h"

void test_crystal_plasticity_single_crystal(Matrix(double) *Fnp1, Matrix(double) *hFnp1, char *filename)
{
  double lame1 = 75600.0;
  double lame2     = 26100.0;
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
  enum {M,pFnp1,xFn,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_redim(double, F2[a],DIM_3,DIM_3);
    Matrix_eye(F2[a],DIM_3);
  } 
  
  double g_n   = mat_p.g0;
  double g_np1 = mat_p.g0;
  double dt     = 0.1;
  double lambda = 0.0;

  staggered_Newton_Rapson_generalized(F2[pFnp1].m_pdata,F2[M].m_pdata, &g_np1, &lambda, 
                                      F2[xFn].m_pdata, F2[xFn].m_pdata,Fnp1->m_pdata, F2[xFn].m_pdata, hFnp1->m_pdata,
                                      g_n, dt, &mat, &elast, &solver_info);
    
  Matrix_print_name(F2[pFnp1], "pFnp1");

  FILE *fp = fopen(filename, "w");  
  fprintf(fp, "%e %e %e\n%e %e %e\n%e %e %e", Mat_v(F2[pFnp1],1,1), Mat_v(F2[pFnp1],1,2), Mat_v(F2[pFnp1],1,3),
                                              Mat_v(F2[pFnp1],2,1), Mat_v(F2[pFnp1],2,2), Mat_v(F2[pFnp1],2,3),
                                              Mat_v(F2[pFnp1],3,1), Mat_v(F2[pFnp1],3,2), Mat_v(F2[pFnp1],3,3));
 
  fclose(fp);  

  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);  

  free(F2);    
  destruct_elasticity(&elast);
  destruct_slip_system(&slip);
}


int main(int argc,char *argv[])
{
  int err = 0;
  
  Matrix(double) Fnp1, hFnp1;
  Matrix_construct_redim(double, Fnp1,  DIM_3, DIM_3);
  Matrix_construct_redim(double, hFnp1, DIM_3, DIM_3);
  Matrix_eye(Fnp1, DIM_3);
  Matrix_eye(hFnp1,DIM_3);
  
  Mat_v(Fnp1,1,2) = Mat_v(Fnp1,1,3) = 0.01;
  test_crystal_plasticity_single_crystal(&Fnp1, &hFnp1, "without_thermal.txt");

  
  Mat_v(hFnp1,1,1) = Mat_v(hFnp1,2,2) = Mat_v(hFnp1,3,3) = 2.0;
   Mat_v(Fnp1,1,1) =  Mat_v(Fnp1,2,2) =  Mat_v(Fnp1,3,3) = 2.0;
   Mat_v(Fnp1,1,2) =  Mat_v(Fnp1,1,3) = 0.02;
  test_crystal_plasticity_single_crystal(&Fnp1, &hFnp1, "with_thermal.txt");
  
  return err;
}
