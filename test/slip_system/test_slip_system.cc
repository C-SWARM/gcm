#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "flowlaw.h"
#include <sys/time.h>

#include <stdio.h>

void test_slip_system(void)
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

  Tensor<2> F = {0.96, 0.00, 0.00,
                 0.00, 0.97, 0.00,
                 0.00, 0.00, 0.97};
  Tensor<2> C;

  MATERIAL_ELASTICITY mat_e;
  set_properties_using_E_and_nu(&mat_e,E,nu);
  print_material_property_elasticity(&mat_e);

  HyperElasticity elast;
  elast.construct_elasticity(&mat_e, true);

  C = F(k,i)*F(k,j);
  elast.update_elasticity(F.data,false);

  SLIP_SYSTEM slip;
  construct_slip_system(&slip,0);

  // create material properties: Plasticity
  MATERIAL_CRYSTAL_PLASTICITY mat_p;
  set_properties_crystal_plasticity(&mat_p,&slip,gamma_dot_0,gamma_dot_s,
                                     m,g0,G0,gs_0,w);
  print_material_property_crystal_plasticity(&mat_p);  // <= this is optional

  // create material plasticity: it needs material properties for elasticity and plasticity
  MATERIAL_CONSTITUTIVE_MODEL mat;
  set_properties_constitutive_model(&mat,&mat_e,&mat_p);

  double *taus         = new double[slip.N_SYS];
  double *gamma_dots   = new double[slip.N_SYS];
  double *dgamma_dtaus = new double[slip.N_SYS];

  double g = 214.4;

  compute_tau_alphas(taus, C.data, elast.S, &slip);
  compute_gamma_dots(gamma_dots, taus, g, &mat_p);
  compute_d_gamma_d_tau(dgamma_dtaus, g, taus, &mat_p);

  printf("eF = [\n%e %e %e \n%e %e %e\n%e %e %e]\n",
          F[0][0], F[0][1], F[0][2],
          F[1][0], F[1][1], F[1][2],
          F[2][0], F[2][1], F[2][2]);

  printf("\n----------------------------------------------\n");
  printf("| alpha | taus | gamma_dots | dgamma_dtaus |\n");
  for(int a=0; a<slip.N_SYS; a++)
    printf("%2d: %e %e %e\n", a+1, taus[a], gamma_dots[a], dgamma_dtaus[a]);

  destruct_slip_system(&slip);

  delete [] taus;
  delete [] gamma_dots;
  delete [] dgamma_dtaus;
}


int main(int argc,char *argv[])
{
  int err = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  test_slip_system();

  gettimeofday(&end, NULL);
  double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0
  + (double)(end.tv_sec - start.tv_sec);
  printf ("Total time: %.4lf s\n", diff);

  return err;
}
