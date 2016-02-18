#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "continuum_damage_model.h"
#include <sys/time.h>

#include <stdio.h>

void test_damage_model(void)
{
  int err = 0;
  
  Matrix(double) F, C, S;
  Matrix_construct_redim(double, F, 3,3);
  Matrix_construct_redim(double, C, 3,3);

  Matrix_eye(F,3);
  
  MATERIAL_ELASTICITY mat_e;
  MATERIAL_CONTINUUM_DAMAGE mat_d;
  
  set_properties_using_E_and_nu(&mat_e,800.0,0.34);
  set_damage_parameters(&mat_d, 8.0, 2.5, 0.15, 100.0, 1.0);
  print_material_property_elasticity(&mat_e);
  
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);
  
  CONTINUUM_DAMAGE damage;

  err += initialize_continuum_damage(&damage,&elast,&mat_d,0.0);  
  
  S.m_row = S.m_col = 3; S.m_pdata = elast.S;  

  double d = 0.001;
  double dt = 1.0;

  for(int a = 0; a<100; a++)
  {
    Mat_v(F,1,1) = 1 + d*a;
    Mat_v(F,2,2) = Mat_v(F,3,3) = 1 - d*a/2.0;

    err += continuum_damage_integration_alg(&damage, dt, F.m_pdata);
    err += update_damage_time_steps(&damage);

    Matrix_AxB(C,1.0,0.0,F,1,F,0);
    update_damaged_elasticity(&damage, dt, F.m_pdata, 1);
    printf("%e, %e, %e\n", d*a, S.m_pdata[0], damage.w);
  }        
  destruct_elasticity(&elast);
  Matrix_cleanup(F);
  Matrix_cleanup(C);
}


int main(int argc,char *argv[])
{
  int err = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  test_damage_model();    
  
  gettimeofday(&end, NULL);
  double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0 
  + (double)(end.tv_sec - start.tv_sec);
  printf ("Total time: %.4lf s\n", diff);
  
  return err;
}
