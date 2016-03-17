#ifndef H__H__SOLVE_SYSTEM__H__H
#define H__H__SOLVE_SYSTEM__H__H

#include "material_properties.h"
#include "hyperelasticity.h"

struct CRYSTAL_PLASTICITY_SOLVER_INFO
{
  int max_itr_stag;     // maximum number of iteration for staggerd NR
  int max_itr_hardening; // maximum number of iteration for hardening NR
  int max_itr_M;         // maximum number of iteration for M NR where (M = pFn_I)
  double tol_hardening;
  double tol_M;
  double computer_zero;
  int max_subdivision;
};

typedef struct CRYSTAL_PLASTICITY_SOLVER_INFO CRYSTAL_PLASTICITY_SOLVER_INFO;

int set_crystal_plasticity_solver_info(CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info, 
                                       int max_itr_stag, int max_itr_hardening, int max_itr_M, 
                                       double tol_hardening, double tol_M, double computer_zero);

int print_crystal_plasticity_solver_info(CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info);     

int staggered_Newton_Rapson(double *pFnp1_out,double *M_out, double *g_out, double *lambda, 
                            double *pFn_in, double *Fn_in, double *Fnp1_in, 
                            double g_n, double dt, 
                            MATERIAL_CONSTITUTIVE_MODEL *mat,
                            ELASTICITY *elasticity, 
                            CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info); 
                            
int Fnp1_Implicit(double *Fnp1_out, double *Fn_in, double *L_in, double dt);

#endif