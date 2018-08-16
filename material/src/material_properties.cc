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

/// print material parameters for the damage model.
///
/// param[in] material property object for the damage model
int print_material_property_damage_model(MATERIAL_CONTINUUM_DAMAGE *mat)
{
  int err = 0;
  printf("-----------------------------------------------------------\n");
  printf("damage model properties\n");
  printf("-----------------------------------------------------------\n");
  printf("Damage scale parameter (p1)     = %e\n", mat->P1);
  printf("Damage shape parameter (p2)     = %e\n", mat->P2);
  printf("Damage energy threshold (Yin)   = %e\n", mat->Y_in);
  printf("Damage viscosity (mu)           = %e\n", mat->mu);  
  printf("Maximum damage                  = %e\n", mat->w_max);
  printf("Split damage parameter alpha_d  = %e\n", mat->alpha_dev);
  printf("Split damage parameter beta_d   = %e\n", mat->beta_dev);
  printf("Split damage parameter alpha_v  = %e\n", mat->alpha_vol);  
  printf("Split damage parameter beta_v   = %e\n", mat->beta_vol);
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
                                          const double M_in,
                                          const double alpha_in,
                                          const double m_in,
                                          const double gamma_dot_0_in,
                                          const double a1_in,
                                          const double a2_in,
                                          const double Lambda1_in,
                                          const double Lambda2_in,
                                          const double c_inf_in,
                                          const double Gamma_in,
                                          const double B_in,
                                          const double pc_b_in,
                                          const double mu_0_in,
                                          const double mu_1_in,
                                          const double p0_in,
                                          const double kappa_in,
                                          const double n_in,
                                          const double g0_in,
                                          const double pc_inf_in)
{                       

  mat->M           = M_in;
  mat->alpha       = alpha_in;
  mat->m           = m_in;
  mat->gamma_dot_0 = gamma_dot_0_in;
  mat->a1          = a1_in;
  mat->a2          = a2_in;
  mat->Lambda1     = Lambda1_in;
  mat->Lambda2     = Lambda2_in;
  mat->c_inf       = c_inf_in;
  mat->Gamma       = Gamma_in;
  mat->B           = B_in;
  mat->pc_b        = pc_b_in;
  mat->mu_0        = mu_0_in;
  mat->mu_1        = mu_1_in;
  mat->p0          = p0_in;
  mat->kappa       = kappa_in;
  mat->n           = n_in;
  mat->g0          = g0_in;
  mat->pc_inf      = pc_inf_in;
}

void print_material_property_poro_visco_plasticity(const MaterialPoroViscoPlasticity *mat)
{
  printf("-----------------------------------------------------------\n");
  printf("Poro-visco-plasticity material properties\n");
  printf("-----------------------------------------------------------\n");
  printf("Yield function parameters      = %e\n", mat->M);      
  printf("  :                            = %e\n", mat->alpha);  
  printf("Flow rule parameters           = %e\n", mat->m);     
  printf("  :                            = %e\n", mat->gamma_dot_0);
  printf("Hardening rule parameters      = %e\n", mat->a1);     
  printf("  :                            = %e\n", mat->a2);     
  printf("  :                            = %e\n", mat->Lambda1);
  printf("  :                            = %e\n", mat->Lambda2);
  printf("Cohesion rule parameters       = %e\n", mat->c_inf);     
  printf("  :                            = %e\n", mat->Gamma);   
  printf("Transition rule parameters     = %e\n", mat->B);       
  printf("  :                            = %e\n", mat->pc_b);     
  printf("Shear modulus parameters       = %e\n", mat->mu_0);      
  printf("  :                            = %e\n", mat->mu_1);      
  printf("Bulk modulus parameters        = %e\n", mat->p0);      
  printf("  :                            = %e\n", mat->kappa);   
  printf("Power law exponent             = %e\n", mat->n);      
  printf("Compaction function parameters = %e\n", mat->g0);     
  printf("  :                            = %e\n", mat->pc_inf);        
}