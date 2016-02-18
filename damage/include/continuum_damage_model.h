#ifndef H__H__CONTINUUM_DAMAGE_MODEL__H__H
#define H__H__CONTINUUM_DAMAGE_MODEL__H__H

#include "material_properties.h"

struct ELASTICITY;
#ifndef TYPE_ELASTICITY
#define TYPE_ELASTICITY
typedef struct ELASTICITY ELASTICITY;
#endif

struct CONTINUUM_DAMAGE;
#ifndef TYPE_CONTINUUM_DAMAGE
#define TYPE_CONTINUUM_DAMAGE
typedef struct CONTINUUM_DAMAGE CONTINUUM_DAMAGE;
#endif

//struct CONTINUUM_DAMAGE_SPLIT;
//#ifndef TYPE_CONTINUUM_DAMAGE_SPLIT
//#define TYPE_CONTINUUM_DAMAGE_SPLIT
//typedef struct CONTINUUM_DAMAGE_SPLIT CONTINUUM_DAMAGE_SPLIT;
//#endif

struct CONTINUUM_DAMAGE
{
  ELASTICITY *elasticity;
  MATERIAL_CONTINUUM_DAMAGE *mat_d;
  double w;
  double X;
  double H;
  double wn;
  double Xn;
  int is_it_damaged;  
};

int initialize_continuum_damage(CONTINUUM_DAMAGE *dam, 
                                ELASTICITY *elast, 
                                MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                double w0);

int continuum_damage_integration_alg(CONTINUUM_DAMAGE *dam,
                                     const double dt,
                                     double *F_in);
                                     
int update_damage_time_steps(CONTINUUM_DAMAGE *dam);

int update_damaged_elasticity(CONTINUUM_DAMAGE *dam, 
                              const double dt,
                              double *F_in, 
                              const int compute_stiffness);                               
//
//struct CONTINUUM_DAMAGE_SPLIT
//{
//  ELASTICITY *elast;
//};

#endif
