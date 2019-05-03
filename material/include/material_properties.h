#ifndef H__H__MATERIAL_PROPERTIES__H__H
#define H__H__MATERIAL_PROPERTIES__H__H

#include<stdio.h>
#include "slip_system.h"

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

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
  double alpha_dev;  // deviatoric energy constibution to deviatoric part of damage
                     // default = 1;
  double beta_dev_p; // volumetirc energy constibution to deviatoric part of damage in tension
                     // default = 1;
  double beta_dev_m; // volumetirc energy constibution to deviatoric part of damage in compression
                     // default = 1;                  
  double alpha_vol;  // deviatoric energy constibution to volume part of damage
                     // default = 1;
  double beta_vol_p; // volumetirc energy constibution to volume part of damage in tension
                     // default = 1;
  double beta_vol_m; // volumetirc energy constibution to volume part of damage in compression
                     // default = 1;                     
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
                                      
int print_material_property_crystal_plasticity(MATERIAL_CRYSTAL_PLASTICITY *mat);                                      

int print_material_property_elasticity(MATERIAL_ELASTICITY *mat);

int set_damage_parameters(MATERIAL_CONTINUUM_DAMAGE *dam, const double P1, 
                                                          const double P2, 
                                                          const double Y_in, 
                                                          const double mu, 
                                                          const double w_max);
                                                          
int set_split_damage_parameters(MATERIAL_CONTINUUM_DAMAGE *dam, const double P1, 
                                                                const double P2, 
                                                                const double Y_in, 
                                                                const double mu, 
                                                                const double w_max,
                                                                const double da,
                                                                const double dbp,
                                                                const double dbm,
                                                                const double va,
                                                                const double vbp,
                                                                const double vbm);

/// print material parameters for the damage model                                                                
int print_material_property_damage_model(MATERIAL_CONTINUUM_DAMAGE *mat);
                                                                
int set_J2_plasticity_parameters(MATERIAL_J2_PLASTICITY *J2P, const double hp,
                                                              const double beta,
                                                              const double k0);
typedef struct
{
  double M,            // Yield function parameters
         alpha,        //   :
         m_d,          // Flow rule parameters
         m_v,
         gamma_dot_0,  //   :
         a1,           // Hardening rule parameters
         a2,           //   :
         Lambda1,      //   :
         Lambda2,      //   :
         c_inf,        // Cohesion rule parameters
         Gamma,        //   :
         B,            // Transition rule parameters
         pc_b,         //   :
         d_m,          //   :
         mu_0,         // Shear modulus parameters
         mu_1,         //   :
         p0,           // Bulk modulus parameters
         kappa,        //   :
         n,            // Power law exponent
         g0,           // Compaction function parameters
         pc_inf,       //   :
         pc_0,         // initial pc
         pJ;           // initial plastic deformation
} MaterialPoroViscoPlasticity;

void set_properties_poro_visco_plasticity(MaterialPoroViscoPlasticity *mat,
                                          const double M_in,
                                          const double alpha_in,
                                          const double m_d_in,
                                          const double m_v_in,                                          
                                          const double gamma_dot_0_in,
                                          const double a1_in,
                                          const double a2_in,
                                          const double Lambda1_in,
                                          const double Lambda2_in,
                                          const double c_inf_in,
                                          const double Gamma_in,
                                          const double B_in,
                                          const double pc_b_in,
                                          const double d_m_in,
                                          const double mu_0_in,
                                          const double mu_1_in,
                                          const double p0_in,
                                          const double kappa_in,
                                          const double n_in,
                                          const double g0_in,
                                          const double pc_inf_in);

void print_material_property_poro_visco_plasticity(const MaterialPoroViscoPlasticity *mat);                                          

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

class MATERIAL_CONSTITUTIVE_MODEL
{
  public:
    MATERIAL_ELASTICITY         *mat_e;
    MATERIAL_CRYSTAL_PLASTICITY *mat_p;
    MATERIAL_CONTINUUM_DAMAGE   *mat_d;
    MATERIAL_J2_PLASTICITY      *mat_J2p;
    MaterialPoroViscoPlasticity *mat_pvp;

  MATERIAL_CONSTITUTIVE_MODEL()
  {
    mat_e   = NULL;
    mat_p   = NULL;
    mat_d   = NULL;
    mat_J2p = NULL;
    mat_pvp = NULL;
  };
};

int set_properties_constitutive_model(MATERIAL_CONSTITUTIVE_MODEL *mat,
                                      MATERIAL_ELASTICITY *mat_e,
                                      MATERIAL_CRYSTAL_PLASTICITY *mat_p);
   
#endif 
