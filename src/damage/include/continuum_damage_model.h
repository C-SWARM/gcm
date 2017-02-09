#ifndef H__H__CONTINUUM_DAMAGE_MODEL__H__H
#define H__H__CONTINUUM_DAMAGE_MODEL__H__H

#include "material_properties.h"

struct ELASTICITY;
#ifndef TYPE_ELASTICITY
#define TYPE_ELASTICITY
typedef struct ELASTICITY ELASTICITY;
#endif

struct CONTINUUM_DAMAGE_SPLIT;
#ifndef TYPE_CONTINUUM_DAMAGE_SPLIT
#define TYPE_CONTINUUM_DAMAGE_SPLIT
typedef struct CONTINUUM_DAMAGE_SPLIT CONTINUUM_DAMAGE_SPLIT;
#endif

int continuum_damage_integration_alg(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                     ELASTICITY *elast,
                                     double *w,
                                     double *X,
                                     double *H,
                                     int *is_it_damaged,
                                     double wn,
                                     double Xn,
                                     const double dt,
                                     double *F_in);

int continuum_damage_integration_alg_public(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                            ELASTICITY *elast,
                                            double *w,
                                            double *X,
                                            double *H,
                                            int *is_it_damaged,
                                            double wn,
                                            double Xn,
                                            const double dt,
                                            double *F_in,
                                            double Y);

int continuum_damage_split_integration_alg(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                           ELASTICITY *elast,
                                           double *dw,
                                           double *vw,
                                           double *dX,
                                           double *vX,
                                           double *dH,
                                           double *vH,
                                           int *is_it_damaged_d,
                                           int *is_it_damaged_v,                                     
                                           double dwn,
                                           double vwn,
                                           double dXn,
                                           double vXn,
                                           const double dt,
                                           double *F_in);
                                           
int continuum_damage_split_integration_alg_public(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                                  ELASTICITY *elast,
                                                  double *dw,
                                                  double *vw,
                                                  double *dX,
                                                  double *vX,
                                                  double *dH,
                                                  double *vH,
                                                  int *is_it_damaged_d,
                                                  int *is_it_damaged_v,                                     
                                                  double dwn,
                                                  double vwn,
                                                  double dXn,
                                                  double vXn,
                                                  const double dt,
                                                  double *F_in,
                                                  double W,
                                                  double U);                                                                                      
                                     
int update_damaged_elasticity(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                              ELASTICITY *elast,
                              double w,
                              int is_it_damaged,
                              double H,
                              const double dt,
                              double *F_in, 
                              const int compute_stiffness);

int update_damaged_elasticity_split(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                    ELASTICITY *elasticity,
                                    double dw,
                                    double vw,
                                    double dH,
                                    double vH,
                                    int is_it_damaged_d,
                                    int is_it_damaged_v,
                                    const double dt,
                                    double *F_in, 
                                    const int compute_stiffness);

int apply_damage_on_stress(double *S, double *S0, double w);
int apply_damage_on_stiffness(double *L_out, double *S0_in, double *L_in, 
                              double w, int is_it_damaged, double H,                               
                              double dt, double mu);
                              
int apply_split_damage_on_stress(double *S, double *dS0, double *vS0, double dw, double vw);
int apply_split_damage_on_stiffness(double *L_out, double *dS0_in, double *vS0_in,
                                    double *dL_in, double *vL_in, 
                                    double dw, double vw, int is_it_damaged_d, int is_it_damaged_v,
                                    double dH, double vH, double dt, double mu);                                   
#endif
