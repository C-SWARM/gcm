#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "J2_plasticity.h"
#include <sys/time.h>

#include <stdio.h>

void test_elasticity(void)
{
  Matrix(double) F, C, S, sp, sp_n, Fn;
  Matrix_construct_init(double, F,    3,3,0.0);
  Matrix_construct_init(double, C,    3,3,0.0);
  Matrix_construct_init(double, sp,   3,3,0.0);
  Matrix_construct_init(double, sp_n, 3,3,0.0);
  Matrix_construct_init(double, Fn,   3,3,0.0);  

  Matrix_eye(F,3);
  Matrix_eye(Fn,3);   
  
  MATERIAL_ELASTICITY mat_e;
  set_properties_using_E_and_nu(&mat_e,70.0e+3,0.25);
  print_material_property_elasticity(&mat_e);
  
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);
  S.m_row = S.m_col = 3; S.m_pdata = elast.S;
  
  double hp = 300.0;
  double beta = 0.0;
  double k0 = 0.5;
  
  double ep = 0.0;
  double ep_n = 0.0;
  double gamma = 0.0;
  double gamma_n = 0.0;
  
  MATERIAL_J2_PLASTICITY J2P;
  set_J2_plasticity_parameters(&J2P,hp,beta,k0);

  double d = 0.001;
  for(int a = 0; a<1000; a++)
  {
    Mat_v(F,1,1) = 1.0 + d*a;
    Mat_v(F,2,2) = Mat_v(F,3,3) = sqrt(1.0/(1.0+d*a));
    Matrix_AxB(C,1.0,0.0,F,1,F,0);
    elast.update_elasticity(&elast,F.m_pdata, 0);

    J2_plasticity_integration_alg(sp.m_pdata,&ep,&gamma,
                                  F.m_pdata,Fn.m_pdata,sp_n.m_pdata,ep_n,
                                  &J2P,&mat_e);


    J2_plasticity_update_elasticity(&J2P,&elast,
                                    F.m_pdata,Fn.m_pdata,sp.m_pdata,sp_n.m_pdata,gamma,0);
                                                                      
    Matrix_AeqB(Fn,1.0,F);
    Matrix_AeqB(sp_n,1.0,sp);
    ep_n = ep;
    gamma_n = gamma;
    
    printf("%e, %e, %e, %e, %e\n", d*a, S.m_pdata[0],gamma, ep, sp.m_pdata[0]);
  }        
  destruct_elasticity(&elast);
  Matrix_cleanup(F);
  Matrix_cleanup(C);
  Matrix_cleanup(sp);
  Matrix_cleanup(sp_n);
  Matrix_cleanup(Fn);
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
