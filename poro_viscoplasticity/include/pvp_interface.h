
#ifndef H_PVP_INTERFACE_H
#define H_PVP_INTERFACE_H

#include "KMS-IJSS2017.h"



int pvp_intf_perform_integration_algorithm2(double *Fnp1,
                                           double *Fn,
                                           double *pFnp1,
                                           double *pFn,
                                           double *pc_np1,
                                           double pc_n,
                                           KMS_IJSS2017_Parameters *mat_pvp,
                                           double dt);

int pvp_intf_perform_integration_algorithm(double *Fnp1,
                                           double *Fn,
                                           double *pFnp1,
                                           double *pFn,
                                           double *pc_np1,
                                           double pc_n,
                                           KMS_IJSS2017_Parameters *mat_pvp,
                                           double dt);

double pvp_intf_hardening_law(double pc,
                              KMS_IJSS2017_Parameters *mat_pvp);

/// compute compaction pressure for given plastic deformation 
double pvp_intf_compute_pc(double pJ, double pc,
                           KMS_IJSS2017_Parameters *mat_pvp);
                           
int pvp_intf_update_elasticity(double *eF,
                               double pc,
                               double *S_in,
                               double *L_in,
                               KMS_IJSS2017_Parameters *mat_pvp,
                               int compute_stiffness);

int pvp_intf_compute_dMdF(double *dMdF_in,
                          double *Fnp1,
                          double *Fn,
                          double *pFnp1,
                          double *pFn,
                          double pc_np1,
                          double pc_n,
                          KMS_IJSS2017_Parameters *mat_pvp,
                          double dt);

#endif

