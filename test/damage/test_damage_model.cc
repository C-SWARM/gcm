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
  Matrix_construct_init(double, C, 3,3,0.0);

  Matrix_eye(F,3);
  
  MATERIAL_ELASTICITY mat_e;
  MATERIAL_CONTINUUM_DAMAGE mat_d;
  
  set_properties_using_E_and_nu(&mat_e,800.0,0.34);
  set_damage_parameters(&mat_d, 8.0, 2.5, 0.15, 100.0, 1.0);
  print_material_property_elasticity(&mat_e);
  
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);
    
  S.m_row = S.m_col = 3; S.m_pdata = elast.S;  

  double d = 0.001;
  double dt = 1.0;
  
  double w, X, H, wn, Xn;
  w = wn = X = Xn = H = 0.0;
  int is_it_damaged = 0;

  FILE *out = fopen("stress.txt", "w");
  for(int a = 0; a<100; a++)
  {
    Mat_v(F,1,1) = 1 + d*a;
    Mat_v(F,2,2) = Mat_v(F,3,3) = 1 - d*a/2.0;

    err += continuum_damage_integration_alg(&mat_d,&elast,
                                            &w,&X,&H,&is_it_damaged,
                                            wn,Xn,dt,F.m_pdata);
    wn = w;
    Xn = X;
    is_it_damaged = 0;

    Matrix_AxB(C,1.0,0.0,F,1,F,0);
    err += update_damaged_elasticity(&mat_d,&elast,w,is_it_damaged,H,
                                     dt,F.m_pdata,1);

    fprintf(out,"%e %e %e\n", d*a, S.m_pdata[0], w);
  }
  fclose(out);        
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
