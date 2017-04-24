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

// class KMS_IJSS2017_Parameters
// *****************************


// Constructors
// ------------

//! This constructor reads and assigns the material parameters
//! It also assigns the smoothMacaulay bool variable, based on which the functions c and d will
//! use a smoothed Macauley brackets or not.
KMS_IJSS2017_Parameters::KMS_IJSS2017_Parameters(double M, double alpha, double m, double gamma_0, double a1, double a2, double L1, double L2,
                                                 double cinf, double Gamma, double B, double pcb, double mu0, double mu1, double p0, double kappa,
                                                 double n, double g0, double pcinf, bool smM)
{
  this->set_parameters(M, alpha, m, gamma_0, a1, a2, L1, L2,
                       cinf, Gamma, B, pcb, mu0, mu1, p0, kappa,
                       n, g0, pcinf, smM);
};


//! This constructor reads and assigns the material parameters
//! It also assigns the smoothMacaulay bool variable, based on which the functions c and d will
//! use a smoothed Macauley brackets or not.
void KMS_IJSS2017_Parameters::set_parameters(double M, double alpha, double m, double gamma_0, double a1, double a2, double L1, double L2,
                                             double cinf, double Gamma, double B, double pcb, double mu0, double mu1, double p0, double kappa,
                                             double n, double g0, double pcinf, bool smM)
{
  
  // assigning the material parameters
  yf_M = M;
  yf_alpha = alpha;
  flr_m = m;
  flr_gamma0=gamma_0;
  hr_a1 = a1;
  hr_a2 = a2;
  hr_Lambda1 = L1;
  hr_Lambda2 = L2;
  c_inf = cinf;
  c_Gamma = Gamma;
  d_B = B;
  d_pcb = pcb;
  mu_0 = mu0;
  mu_1 = mu1;
  K_p0 = p0;
  K_kappa = kappa;
  pl_n = n;
  cf_g0 = g0;
  cf_pcinf = pcinf;
  
  // smooth Macauley brackets
  smMbrackets = smM;
};


void KMS_IJSS2017_Parameters::Checks( bool Verbose )
//! This method performs some checks and prints warnings, foreseeing issues in the
//! KMS_IJSS2017 model integration
{
  if ( Verbose && ( hr_a1 * hr_Lambda1 <=0 || hr_a2 * hr_Lambda2 <=0 ) )
  {
    std::cout << " WARNING: KMS_IJSS2017_Explicit<dim>::FindpcFromJp( ) - The sign of the first derivative is not defined.";
    std::cout << " Code did not abort but outcomes might be wrong. \n";
  };
};




// IO - methods
// ------------

void KMS_IJSS2017_Parameters::AsAString( std::string& str)
//! This method prints the model features as a string
{
  
  str += "\n";
  str += "  KMS_IJSS2017 model material parameters:\n";
  str += "   Yield function parameters:\n";
  str += "    M = " + std::to_string(yf_M);
  str += ", alpha = " + std::to_string(yf_alpha);
  str += "\n   Flow rule parameters:\n";
  str += "    m = " + std::to_string(flr_m);
  str += ", gamma0 = " + std::to_string(flr_gamma0);
  str += "\n   Hardening rule parameters:\n";
  str += "    a1 = " + std::to_string(hr_a1);
  str += ", a2 = " + std::to_string(hr_a2);
  str += ", Lambda1 = " + std::to_string(hr_Lambda1);
  str += ", Lambda2 = " + std::to_string(hr_Lambda2);
  str += "\n   Cohesion rule parameters:\n";
  str += "    c_inf = " + std::to_string(c_inf);
  str += ", Gamma = " + std::to_string(c_Gamma);
  str += "\n   Transition rule parameters:\n";
  str += "    B = " + std::to_string(d_B);
  str += ", pcb = " + std::to_string(d_pcb);
  str += "\n   Shear rule parameters:\n";
  str += "    mu_0 = " + std::to_string(mu_0);
  str += ", mu_1 = " + std::to_string(mu_1);
  str += "\n   Bulk rule parameters:\n";
  str += "    p_0 = " + std::to_string(K_p0);
  str += ", kappa = " + std::to_string(K_kappa);
  str += "\n   Power law exponent:\n";
  str += "    n = " + std::to_string(pl_n);
  str += "\n   Compaction function parameters:\n";
  str += "    g0 = " + std::to_string(cf_g0);
  str += ", pcinf = " + std::to_string(cf_pcinf);
  
  str += "\n";
};                                                              