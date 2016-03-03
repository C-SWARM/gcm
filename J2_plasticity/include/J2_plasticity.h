#ifndef H__H__J2_PLASTICITY__H__H
#define H__H__J2_PLASTICITY__H__H

struct MATERIAL_ELASTICITY;
#ifndef TYPE_MATERIAL_ELASTICITY
#define TYPE_MATERIAL_ELASTICITY
typedef struct MATERIAL_ELASTICITY MATERIAL_ELASTICITY;
#endif

struct MATERIAL_J2_PLASTICITY;
#ifndef TYPE_MATERIAL_J2_PLASTICITY
#define TYPE_MATERIAL_J2_PLASTICITY
typedef struct MATERIAL_J2_PLASTICITY MATERIAL_J2_PLASTICITY;
#endif

struct ELASTICITY;
#ifndef TYPE_ELASTICITY
#define TYPE_ELASTICITY
typedef struct ELASTICITY ELASTICITY;
#endif

int J2_plasticity_integration_alg(double *sp_out,
                                  double *ep_out,
                                  double *gamma_out,
                                  double *F_in, 
                                  double *Fn_in,
                                  double *sp_n_in,
                                  double ep_n,
                                  MATERIAL_J2_PLASTICITY *J2P,
                                  MATERIAL_ELASTICITY *mat_e);
                                  

int J2_plasticity_update_elasticity(MATERIAL_J2_PLASTICITY *J2P,
                                    ELASTICITY *elast,
                                    double *F_in,
                                    double *Fn_in,
                                    double *sp_in,
                                    double *sp_n_in,                                    
                                    double gamma,
                                    const int compute_stiffness);
#endif