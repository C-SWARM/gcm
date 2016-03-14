#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "J2_plasticity.h"
#include <sys/time.h>

#include <stdio.h>

void test_elasticity(void)
{
  Matrix(double) F, FT, C, S, E, e, sp, sp_n, Fn, sigma;
  Matrix_construct_init(double, F,    3,3,0.0);
  Matrix_construct_init(double, FT,   3,3,0.0);  
  Matrix_construct_init(double, C,    3,3,0.0);
  Matrix_construct_init(double, E,    3,3,0.0);
  Matrix_construct_init(double, e,    3,3,0.0);
  Matrix_construct_init(double, sp,   3,3,0.0);
  Matrix_construct_init(double, sp_n, 3,3,0.0);
  Matrix_construct_init(double, Fn,   3,3,0.0);
  Matrix_construct_init(double, sigma,   3,3,0.0);  
  
  Matrix_eye(F,3);
  Matrix_eye(Fn,3);   
  
  MATERIAL_ELASTICITY mat_e;
  set_properties_using_E_and_nu(&mat_e,9.95e+3,0.25);
  print_material_property_elasticity(&mat_e);
  
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);
  S.m_row = S.m_col = 3; S.m_pdata = elast.S;
  
  double hp = 0.2;
  double beta = 0.0;
  double k0 = 200.0;
  
  double ep = 0.0;
  double ep_n = 0.0;
  double gamma = 0.0;
  double gamma_n = 0.0;
  
  MATERIAL_J2_PLASTICITY J2P;
  set_J2_plasticity_parameters(&J2P,hp,beta,k0);

  FILE *out = fopen("out.txt", "w");
  double d = 0.001;
  for(int a = 0; a<220; a++)
  {
    Mat_v(F,1,1) = 1.0 - d*a;
    Mat_v(F,2,2) = Mat_v(F,3,3) = 1.0/sqrt(1.0 - 0.5*d*a);
    Matrix_AxB(C,1.0,0.0,F,1,F,0);
    Matrix_AeqB(E,0.5,C);
    Mat_v(E,1,1) -= 0.5;
    Mat_v(E,2,2) -= 0.5;
    Mat_v(E,3,3) -= 0.5;
    
    Matrix_AeqBT(FT,1.0,F);
    Matrix_Tns2_AxBxC(e,1.0,0.0,F,E,FT);
        
    elast.update_elasticity(&elast,F.m_pdata, 0);

    J2_plasticity_integration_alg(sp.m_pdata,&ep,&gamma,
                                  F.m_pdata,Fn.m_pdata,sp_n.m_pdata,ep_n,
                                  &J2P,&mat_e);

    
    J2_plasticity_update_elasticity(&J2P,&elast,
                                    F.m_pdata,Fn.m_pdata,sp.m_pdata,sp_n.m_pdata,gamma,0);
    
    double sigma_eff = 0.0;
    elast.compute_Cauchy_eff(&elast,&sigma_eff,F.m_pdata);
    elast.compute_Cauchy(&elast,sigma.m_pdata,F.m_pdata);
                                                                          
    Matrix_AeqB(Fn,1.0,F);
    Matrix_AeqB(sp_n,1.0,sp);
    ep_n = ep;
    gamma_n = gamma;
    
    double et = log(1.0-d*a);
    fprintf(out, "%e, %e, %e, %e, %e, %e, %e\n", et, 
                                                 Mat_v(sigma,1,1),
                                                 (elast.S[0]),
                                                 sigma_eff, 
                                                 gamma, 
                                                 ep, 
                                                 sp.m_pdata[0]);
  }
  fclose(out);
          
  destruct_elasticity(&elast);
  Matrix_cleanup(F);
  Matrix_cleanup(FT);  
  Matrix_cleanup(C);
  Matrix_cleanup(E);
  Matrix_cleanup(e);
  Matrix_cleanup(sp);
  Matrix_cleanup(sp_n);
  Matrix_cleanup(Fn);
  Matrix_cleanup(sigma);  
}


int main(int argc,char *argv[])
{
  int err = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  test_elasticity();    
  
  gettimeofday(&end, NULL);
  double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0 
  + (double)(end.tv_sec - start.tv_sec);
  printf ("Total time: %.4lf s\n", diff);
  
  return err;
}
