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
  dam->P1      = P1;
  dam->P2      = P2;
  dam->Y_in    = Y_in;
  dam->mu      = mu;
  dam->w_max   = w_max;
  dam->alpha_dev = 1.0;
  dam->beta_dev  = 0.0;
  dam->alpha_vol = 0.0;
  dam->beta_vol  = 1.0;
  return err;
}

int set_split_damage_parameters(MATERIAL_CONTINUUM_DAMAGE *dam, const double P1, 
                                                                const double P2, 
                                                                const double Y_in, 
                                                                const double mu, 
                                                                const double w_max,
                                                                const double da,
                                                                const double db,
                                                                const double va,
                                                                const double vb)
{
  int err = 0;
  dam->P1      = P1;
  dam->P2      = P2;
  dam->Y_in    = Y_in;
  dam->mu      = mu;
  dam->w_max   = w_max;
  dam->alpha_dev = da;
  dam->beta_dev  = db;
  dam->alpha_vol = va;
  dam->beta_vol  = vb;
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

void set_properties_poro_visco_plasticity(MaterialPoroViscoPlasticity *mat,
                                          const double yf_M_in,
                                          const double yf_alpha_in,
                                          const double flr_m_in,
                                          const double flr_gamma_dot_0_in,
                                          const double hr_a1_in,
                                          const double hr_a2_in,
                                          const double hr_Lambda1_in,
                                          const double hr_Lambda2_in,
                                          const double c_inf_in,
                                          const double c_Gamma_in,
                                          const double d_B_in,
                                          const double d_pcb_in,
                                          const double mu_0_in,
                                          const double mu_1_in,
                                          const double K_p0_in,
                                          const double K_kappa_in,
                                          const double pl_n_in,
                                          const double cf_g0_in,
                                          const double cf_pcinf_in)
{                       

  mat->yf_M            = yf_M_in;
  mat->yf_alpha        = yf_alpha_in;
  mat->flr_m           = flr_m_in;
  mat->flr_gamma_dot_0 = flr_gamma_dot_0_in;
  mat->hr_a1           = hr_a1_in;
  mat->hr_a2           = hr_a2_in;
  mat->hr_Lambda1      = hr_Lambda1_in;
  mat->hr_Lambda2      = hr_Lambda2_in;
  mat->c_inf           = c_inf_in;
  mat->c_Gamma         = c_Gamma_in;
  mat->d_B             = d_B_in;
  mat->d_pcb           = d_pcb_in;
  mat->mu_0            = mu_0_in;
  mat->mu_1            = mu_1_in;
  mat->K_p0            = K_p0_in;
  mat->K_kappa         = K_kappa_in;
  mat->pl_n            = pl_n_in;
  mat->cf_g0           = cf_g0_in;
  mat->cf_pcinf        = cf_pcinf_in;
}

void print_material_property_poro_visco_plasticity(const MaterialPoroViscoPlasticity *mat)
{
  printf("-----------------------------------------------------------\n");
  printf("Poro-visco-plasticity material properties\n");
  printf("-----------------------------------------------------------\n");
  printf("Yield function parameters      = %e\n", mat->yf_M);      
  printf("  :                            = %e\n", mat->yf_alpha);  
  printf("Flow rule parameters           = %e\n", mat->flr_m);     
  printf("  :                            = %e\n", mat->flr_gamma_dot_0);
  printf("Hardening rule parameters      = %e\n", mat->hr_a1);     
  printf("  :                            = %e\n", mat->hr_a2);     
  printf("  :                            = %e\n", mat->hr_Lambda1);
  printf("  :                            = %e\n", mat->hr_Lambda2);
  printf("Cohesion rule parameters       = %e\n", mat->c_inf);     
  printf("  :                            = %e\n", mat->c_Gamma);   
  printf("Transition rule parameters     = %e\n", mat->d_B);       
  printf("  :                            = %e\n", mat->d_pcb);     
  printf("Shear modulus parameters       = %e\n", mat->mu_0);      
  printf("  :                            = %e\n", mat->mu_1);      
  printf("Bulk modulus parameters        = %e\n", mat->K_p0);      
  printf("  :                            = %e\n", mat->K_kappa);   
  printf("Power law exponent             = %e\n", mat->pl_n);      
  printf("Compaction function parameters = %e\n", mat->cf_g0);     
  printf("  :                            = %e\n", mat->cf_pcinf);        
}