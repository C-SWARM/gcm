#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "crystal_plasticity_integration.h"
#include "flowlaw.h"
#include <sys/time.h>

void test_crystal_plasticity_single_crystal(void)
{
  double lame1 = 75600.0;
  double lame2     = 26100.0;
  double E = 207.0e+3;
  double nu = 0.31;
  
  double gamma_dot_0 = 1.0;
  double gamma_dot_s = 50.0e+9;
  double m           = 0.06;  
  double g0          = 45.0;
  double G0          = 500.0;
  double gs_0        = 155;
  double w           = 5.0e-4;

  int max_itr_stag      = 100;
  int max_itr_hardening = 5;
  int max_itr_M         = 100;
  double tol_hardening  = 1.0e-5;
  double tol_M          = 1.0e-5;
  double computer_zero  = 1.0e-15;

  // create material properties: Elasticity
  MATERIAL_ELASTICITY mat_e;  
  set_properties_using_E_and_nu(&mat_e,E,nu);
  mat_e.m01 =   mat_e.m10 = 18.254e+3;
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
  enum {M,MI,pFn,pFnp1,pFnp1_I,eFnp1,Fn,Fnp1,L,sigma,PK2dev,sigma_dev,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_init(double, F2[a],DIM_3,DIM_3,0.0);
    Matrix_eye(F2[a],DIM_3);
  } 
  
  double g_n,g_np1;
  g_n = g_np1 = mat_p.g0;
  
  double dt = 0.001;
  
  double d = 1.0;
  // set velocity gradient  
  Mat_v(F2[L],1,1) = -d;
  Mat_v(F2[L],2,2) = Mat_v(F2[L],3,3) = d/2;  
  
  // start integration  
  Matrix(double) PK2;
  PK2.m_row = PK2.m_col = DIM_3; PK2.m_pdata = elast.S;
  
  FILE *fp = fopen("single_crystal_results.txt", "w");

F2[Fnp1].m_pdata[0] = 9.999924e-01; F2[Fnp1].m_pdata[1] = -1.819720e-06; F2[Fnp1].m_pdata[2] = -1.523077e-07; 
F2[Fnp1].m_pdata[3] =-2.542806e-06; F2[Fnp1].m_pdata[4] = 9.999904e-01; F2[Fnp1].m_pdata[5] = -3.102639e-07; 
F2[Fnp1].m_pdata[6] =3.789849e-06; F2[Fnp1].m_pdata[7] = 2.592934e-06; F2[Fnp1].m_pdata[8] = 1.000026e+00;
  
  for(int a = 1; a<=1; a++)
  {
    double lambda = 0.0;
    double t = a*dt;
    
    // compute total deformation gradient using velocity gradient
    //Fnp1_Implicit(F2[Fnp1].m_pdata, F2[Fn].m_pdata, F2[L].m_pdata, dt); 
    
    staggered_Newton_Rapson(F2[pFnp1].m_pdata,F2[M].m_pdata, &g_np1, &lambda, 
                            F2[pFn].m_pdata, F2[Fn].m_pdata,F2[Fnp1].m_pdata, 
                            g_n, dt, &mat, &elast, &solver_info);
    Matrix_AeqB(F2[pFn],1.0,F2[pFnp1]);
    Matrix_AeqB(F2[Fn],1.0,F2[Fnp1]);  
    Matrix_inv(F2[pFnp1],F2[pFnp1_I]);
    Matrix_AxB(F2[eFnp1],1.0,0.0,F2[Fnp1],0,F2[pFnp1_I],0);    
    g_n = g_np1;
    
    
    // print result at time t
    double sigma_eff;
    double PK2_eff;
    double det_pF;
    Matrix_det(F2[pFnp1], det_pF);
    
    elast.update_elasticity(&elast,F2[eFnp1].m_pdata,0);
    elast.compute_PK2_eff(&elast,&PK2_eff);
    elast.compute_Cauchy_eff(&elast,&sigma_eff,F2[eFnp1].m_pdata);   

    fprintf(fp, "%e %e %e %e %e %e\n",t,sigma_eff,PK2_eff, g_np1, det_pF, Mat_v(PK2,1,1));
  }    
  
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
  struct timeval start, end;
  gettimeofday(&start, NULL);

  test_crystal_plasticity_single_crystal();
  
  gettimeofday(&end, NULL);
  double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0 
  + (double)(end.tv_sec - start.tv_sec);
  printf ("Total time: %.4lf s\n", diff);
  
  return err;
}
