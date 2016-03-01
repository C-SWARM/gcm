#ifndef H__H__MATERIAL_PROPERTIES__H__H
#define H__H__MATERIAL_PROPERTIES__H__H

#include "slip_system.h"

struct MATERIAL_ELASTICITY
{
  double E;       // Young's modulus
  double nu;      // Poission's ratio
  double G;       // Shear modulus
  double lambda;  // 1st Lame constant
  double mu;      // 2nd Lame constant (shear modulus)  
  double m01;     // Mooney Rivlin parameter
  double m10;     // Mooney Rivlin parameter
  double kappa;   // bulk modulus
  int volPotFlag; // potential flag (volume)
  int devPotFlag; // potential flag (deviatoric)
};

struct MATERIAL_CRYSTAL_PLASTICITY
{
  SLIP_SYSTEM *slip;
  double gamma_dot_0; // reference rate
  double gamma_dot_s; //
  double m;           // rate sensitivity  
  double g0;          // initial resolved shear strength
  double G0;          // hardening rate
  double gs_0;        
  double w;
};

struct MATERIAL_CONTINUUM_DAMAGE
{
  double P1;      // Weibull function parameters: scale
  double P2;      // Weibull function parameters: shape
  double Y_in;    // Weibull function parameters: initial damage energy threshold
  double mu;      // damage viscosity
  double w_max;   // maximum damage
};

struct MATERIAL_J2_PLASTICITY
{
  double hp;      // 
  double beta;    // 
  double k0;      // 
};

typedef struct MATERIAL_ELASTICITY MATERIAL_ELASTICITY;
typedef struct MATERIAL_CRYSTAL_PLASTICITY MATERIAL_CRYSTAL_PLASTICITY;
typedef struct MATERIAL_CONTINUUM_DAMAGE MATERIAL_CONTINUUM_DAMAGE;
typedef struct MATERIAL_J2_PLASTICITY MATERIAL_J2_PLASTICITY;

struct MATERIAL_CONSTITUTIVE_MODEL
{
  MATERIAL_ELASTICITY *mat_e;
  MATERIAL_CRYSTAL_PLASTICITY *mat_p;
};
typedef struct MATERIAL_CONSTITUTIVE_MODEL MATERIAL_CONSTITUTIVE_MODEL;

/**
 * set material properties for elasticity (Mooney_Rivlin and Doll_Schweizerhof_7)
 * \param[in] E and nu
 * \param[out] mat
 * \return non-zero on internal error.
 */
int set_properties_using_E_and_nu(MATERIAL_ELASTICITY *mat, double E, double nu);

/**
 * set material properties for elasticity (Mooney_Rivlin and Doll_Schweizerhof_7)
 * \param[in] lambda and mu
 * \param[out] mat
 * \return non-zero on internal error.
 */
int set_properties_using_Lame_constants(MATERIAL_ELASTICITY *mat, double lambda, double mu);

/**
 * set material properties for elasticity for Mooney_Rivlin
 * \param[in] E, nu, m01, m10, and volPotFlag
 * \param[out] mat
 * \return non-zero on internal error.
 */
int set_properties_Mooney_Rivlin(MATERIAL_ELASTICITY *mat, double E, double nu, double m01, double m10, int volPotFlag);

/**
 * set material properties for Linear elasticity
 * \param[in] E, nu, and volPotFlag
 * \param[out] mat
 * \return non-zero on internal error.
 */
int set_properties_simple_linear(MATERIAL_ELASTICITY *mat, double E, double nu, int volPotFlag);
int set_properties_crystal_plasticity(MATERIAL_CRYSTAL_PLASTICITY *mat, 
                                      SLIP_SYSTEM *slip, double gamma_dot_0, double gamma_dot_s, 
                                      double m, double g0, double G0, double gs_0, double w);
                                      
int set_properties_constitutive_model(MATERIAL_CONSTITUTIVE_MODEL *mat,
                                      MATERIAL_ELASTICITY *mat_e,
                                      MATERIAL_CRYSTAL_PLASTICITY *mat_p);

int print_material_property_crystal_plasticity(MATERIAL_CRYSTAL_PLASTICITY *mat);                                      

int print_material_property_elasticity(MATERIAL_ELASTICITY *mat);

int set_damage_parameters(MATERIAL_CONTINUUM_DAMAGE *dam, const double P1, 
                                                          const double P2, 
                                                          const double Y_in, 
                                                          const double mu, 
                                                          const double w_max);

int set_J2_plasticity_parameters(MATERIAL_J2_PLASTICITY *J2P, const double hp,
                                                              const double beta,
                                                              const double k0);
#endif 