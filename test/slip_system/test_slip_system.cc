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
  double m_m[9] = {0.96, 0.00, 0.00, 
                   0.00, 0.97, 0.00, 
                   0.00, 0.00, 0.97};
  
  double gamma_dot_0 = 1.0;
  double gamma_dot_s = 50.0e+9;
  double m           = 0.05;  
  double g0          = 210.0;
  double G0          = 200.0;
  double gs_0        = 330.0;
  double w           = 0.005;
    
  Matrix(double) F, C;  
  Matrix_construct_init(double, C, 3,3,0.0);  
  
  F.m_row = F.m_col = 3;
  F.m_pdata = m_m;
  
  MATERIAL_ELASTICITY mat_e;
  set_properties_using_E_and_nu(&mat_e,E,nu);
  print_material_property_elasticity(&mat_e);
  
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);

  Matrix_AxB(C,1.0,0.0,F,1,F,0);
  elast.update_elasticity(&elast,F.m_pdata,0);

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
    
  
  Matrix(double) taus, gamma_dots, dgamma_dtaus;
  Matrix_construct_redim(double, taus, slip.N_SYS,1);
  Matrix_construct_redim(double, gamma_dots, slip.N_SYS,1);
  Matrix_construct_redim(double, dgamma_dtaus, slip.N_SYS,1);
  
  double g = 214.4;
  
  compute_tau_alphas(taus.m_pdata, C.m_pdata, elast.S, &slip);
  compute_gamma_dots(gamma_dots.m_pdata, taus.m_pdata, g, &mat_p);  
  compute_d_gamma_d_tau(dgamma_dtaus.m_pdata, g, taus.m_pdata, &mat_p);
  
  Matrix_print_name(F, "eF");
  
  printf("\n----------------------------------------------\n");
  printf("| alpha | taus | gamma_dots | dgamma_dtaus |\n");
  for(int a=0; a<slip.N_SYS; a++)
    printf("%2d: %e %e %e\n", a+1, taus.m_pdata[a], gamma_dots.m_pdata[a], dgamma_dtaus.m_pdata[a]);
  
  destruct_slip_system(&slip);
  destruct_elasticity(&elast);
  
  Matrix_cleanup(C); 
  Matrix_cleanup(taus);
  Matrix_cleanup(gamma_dots);
  Matrix_cleanup(dgamma_dtaus);          
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
