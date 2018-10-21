#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include <sys/time.h>

#include <stdio.h>

void test_elasticity(void)
{
  Tensor<2> F = {}, C;
  
  MATERIAL_ELASTICITY mat_e;
  set_properties_using_E_and_nu(&mat_e,70.0e+3,0.25);
  print_material_property_elasticity(&mat_e);
  
  GcmElasticity<MATERIAL_ELASTICITY, StrainEnergyDensityFunction, void *> elast(&mat_e, true);
  TensorA<2> S(elast.S);

  double d = 0.01;
  FILE *out = fopen("stress.txt", "w");
  for(int a = 0; a<100; a++)
  {
    F[0][0] = 1.0 + d*a;
    F[1][1] = F[2][2] = 1.0 - d*a/2.0;
    C = F(k,i)*F(k,j);
    elast.update_elasticity(F.data, false);
    fprintf(out, "%e %e\n", d*a, S[0][0]);
  }
  fclose(out);
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
