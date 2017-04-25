#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include <sys/time.h>

#include <stdio.h>

void test_elasticity(void)
{
  Matrix(double) F, C, S;
  Matrix_construct_redim(double, F, 3,3);
  Matrix_construct_init(double, C, 3,3,0.0);

  Matrix_eye(F,3);
  
  MATERIAL_ELASTICITY mat_e;
  set_properties_using_E_and_nu(&mat_e,70.0e+3,0.25);
  print_material_property_elasticity(&mat_e);
  
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);
  S.m_row = S.m_col = 3; S.m_pdata = elast.S;  

  double d = 0.01;
  FILE *out = fopen("stress.txt", "w");
  for(int a = 0; a<100; a++)
  {
    Mat_v(F,1,1) = 1 + d*a;
    Mat_v(F,2,2) = Mat_v(F,3,3) = 1 - d*a/2.0;
    Matrix_AxB(C,1.0,0.0,F,1,F,0);
    elast.update_elasticity(&elast,F.m_pdata, 0);
    fprintf(out, "%e %e\n", d*a, S.m_pdata[0]);
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

  test_elasticity();    
  
  gettimeofday(&end, NULL);
  double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0 
  + (double)(end.tv_sec - start.tv_sec);
  printf ("Total time: %.4lf s\n", diff);
  
  return err;
}
