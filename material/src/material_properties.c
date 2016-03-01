#include "material_properties.h"
#include <stdio.h>

int set_properties_using_E_and_nu(MATERIAL_ELASTICITY *mat, double E, double nu)
{
  int err = 0;
  mat->E      = E;
  mat->nu     = nu;
  mat->lambda = E*nu/(1.0+nu)/(1.0-2.0*nu);
  mat->mu     = E/2.0/(1.0+nu);
  mat->G      = mat->mu;  
  mat->kappa  = E/3.0/(1.0-2.0*nu);
  mat->m01    = 0.0;         // mu = 2(m10 + m01), For consistency with linear elasticity
  mat->m10    = mat->mu/2.0; // Neo-Hookean
  mat->devPotFlag = 1;
  mat->volPotFlag = 2;  
  return err;    
}

int set_properties_using_Lame_constants(MATERIAL_ELASTICITY *mat, double lambda, double mu)
{
  int err = 0;
  double E  = mu*((3.0*lambda+2.0*mu)/(lambda+mu));
  double nu = E/(2.0*mu)-1.0;
  mat->E      = E;
  mat->nu     = nu;
  mat->mu     = mu;
  mat->G      = mat->mu;  
  mat->lambda = lambda;
  mat->kappa  = E/3.0/(1.0-2.0*nu);
  mat->m01    = 0.0;         // mu = 2(m10 + m01), For consistency with linear elasticity
  mat->m10    = mat->mu/2.0; // Neo-Hookean
  mat->devPotFlag = 1;
  mat->volPotFlag = 2;   
  return err;    
}

int set_properties_Mooney_Rivlin(MATERIAL_ELASTICITY *mat, double E, double nu, double m01, double m10, int volPotFlag)
{
  int err = 0;
  mat->E      = E;
  mat->nu     = nu;
  mat->lambda = E*nu/(1.0+nu)/(1.0-2.0*nu);
  mat->mu     = E/2.0/(1.0+nu);
  mat->G      = mat->mu;  
  mat->kappa  = E/3.0/(1.0-2.0*nu);
  mat->m01    = m01;
  mat->m10    = m10;
  mat->devPotFlag = 1;
  mat->volPotFlag = volPotFlag;   
  return err;    
}

int set_properties_simple_linear(MATERIAL_ELASTICITY *mat, double E, double nu, int volPotFlag)
{
  int err = 0;
  mat->E      = E;
  mat->nu     = nu;
  mat->lambda = E*nu/(1.0+nu)/(1.0-2.0*nu);
  mat->mu     = E/2.0/(1.0+nu);
  mat->G      = mat->mu;  
  mat->kappa  = E/3.0/(1.0-2.0*nu);
  mat->m01    = 0.0;         // mu = 2(m10 + m01), For consistency with linear elasticity
  mat->m10    = mat->mu/2.0; // Neo-Hookean
  mat->devPotFlag = 2;
  mat->volPotFlag = volPotFlag;   
  return err;    
}

int set_properties_crystal_plasticity(MATERIAL_CRYSTAL_PLASTICITY *mat, 
                                      SLIP_SYSTEM *slip, double gamma_dot_0, double gamma_dot_s, 
                                      double m, double g0, double G0, double gs_0, double w)
{
  int err = 0;
  mat->slip        = slip;
  mat->gamma_dot_0 = gamma_dot_0;
  mat->gamma_dot_s = gamma_dot_s;
  mat->m           = m;
  mat->g0          = g0;
  mat->G0          = G0;
  mat->gs_0        = gs_0;
  mat->w           = w;
  return err;              
}

int set_properties_constitutive_model(MATERIAL_CONSTITUTIVE_MODEL *mat,
                                      MATERIAL_ELASTICITY *mat_e,
                                      MATERIAL_CRYSTAL_PLASTICITY *mat_p)
{
  int err = 0;
  mat->mat_e = mat_e;
  mat->mat_p = mat_p;
  return err;
}

int print_material_property_crystal_plasticity(MATERIAL_CRYSTAL_PLASTICITY *mat)
{
  int err = 0;
  printf("-----------------------------------------------------------\n");
  printf("crystal plasticity material properties\n");
  printf("-----------------------------------------------------------\n");      
  printf("%s, N_SYS = %d\n",   mat->slip->name, mat->slip->N_SYS);  
  printf("gamma_dot_0 = %e\n", mat->gamma_dot_0);
  printf("gamma_dot_s = %e\n", mat->gamma_dot_s);
  printf("m           = %e\n", mat->m);
  printf("g0          = %e\n", mat->g0);
  printf("G0          = %e\n", mat->G0);
  printf("gs_0        = %e\n", mat->gs_0);
  printf("w           = %e\n", mat->w);
  return err;
}
                                                 
int print_material_property_elasticity(MATERIAL_ELASTICITY *mat)
{
  int err = 0;
  printf("-----------------------------------------------------------\n");
  printf("elasticity material properties\n");
  printf("-----------------------------------------------------------\n");
  printf("Young's modulus (E)                    = %e\n", mat->E);
  printf("Poission ratio (nu)                    = %e\n", mat->nu);
  printf("Shear modulus (G)                      = %e\n", mat->G);
  printf("Lame constant 1 (lambda)               = %e\n", mat->lambda);  
  printf("Lame constant 2 (mu)                   = %e\n", mat->mu);
  printf("Bulk modulus (kappa)                   = %e\n", mat->kappa);
  printf("Mooney_Rivlin parameters (m01)         = %e\n", mat->m01);
  printf("Mooney_Rivlin parameters (m10)         = %e\n", mat->m10);  
  printf("SEDF flag deviatoric part (devPotFlag) = %d\n", mat->devPotFlag);
  printf("SEDF flag volume part (volPotFlag)     = %d\n", mat->volPotFlag);   
  return err;
}

int set_damage_parameters(MATERIAL_CONTINUUM_DAMAGE *dam, const double P1, 
                                                          const double P2, 
                                                          const double Y_in, 
                                                          const double mu, 
                                                          const double w_max)
{
  int err = 0;
  dam->P1    = P1;
  dam->P2    = P2;
  dam->Y_in  = Y_in;
  dam->mu    = mu;
  dam->w_max = w_max;
  return err;
}

int set_J2_plasticity_parameters(MATERIAL_J2_PLASTICITY *J2P, const double hp,
                                                              const double beta,
                                                              const double k0)
{
  int err = 0;
  J2P->hp = hp;
  J2P->beta = beta;
  J2P->k0 = k0;
  return err;
}                                                              