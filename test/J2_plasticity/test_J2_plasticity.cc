#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "J2_plasticity.h"
#include <sys/time.h>

#include <stdio.h>

void test_elasticity(void)
{
  Tensor<2> F, eye, E, e, sp, sp_n, Fn, sigma;
  
  eye = ttl::identity(i,j); 
  F  = eye(i,j);
  Fn = eye(i,j);   
  
  MATERIAL_ELASTICITY mat_e;
  set_properties_using_E_and_nu(&mat_e,9.95e+3,0.25);
  print_material_property_elasticity(&mat_e);
  
  HyperElasticity elast;
  elast.construct_elasticity(&mat_e, true);
  
  double hp = 0.2;
  double beta = 0.0;
  double k0 = 200.0;
  
  double ep = 0.0;
  double ep_n = 0.0;
  double gamma = 0.0;
  
  MATERIAL_J2_PLASTICITY J2P;
  set_J2_plasticity_parameters(&J2P,hp,beta,k0);

  FILE *out = fopen("stress.txt", "w");
  double d = 0.001;
  for(int a = 0; a<220; a++)
  {
    F[0][0] = 1.0 - d*a;
    F[1][1] = F[2][2] = 1.0/sqrt(1.0 - 0.5*d*a);
    
     
    E = 0.5*(F(k,i)*F(k,j) - eye(i,j));

    
    e = F(i,k)*E(k,l)*F(j,l);
            
    elast.update_elasticity(F.data, false);

    J2_plasticity_integration_alg(sp.data,&ep,&gamma,
                                  F.data,Fn.data,sp_n.data,ep_n,
                                  &J2P,&mat_e);

    
    J2_plasticity_update_elasticity(&J2P,&elast,
                                    F.data,Fn.data,sp.data,sp_n.data,gamma,0);
    
    double sigma_eff = 0.0;
    elast.compute_Cauchy_eff(&sigma_eff,F.data);
    elast.compute_Cauchy(sigma.data,F.data);

    Fn   = F(i,j);
    sp_n = sp(i,j);
    
    ep_n = ep;
    
    double et = log(1.0-d*a);
    fprintf(out, "%e %e %e %e %e %e %e\n", et, 
                                           sigma[0][0],
                                           elast.S[0],
                                           sigma_eff, 
                                           gamma, 
                                           ep, 
                                           sp[0][0]);
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
