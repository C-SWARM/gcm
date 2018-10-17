#include"GcmSolverInfo.h"
#include<stdio.h>

void set_gcm_solver_info(GcmSolverInfo *solver_info,
                         const int max_itr_stag, 
                         const int max_itr_hardening, 
                         const int max_itr_M,
                         const double tol_hardening, 
                         const double tol_M,
                         const double computer_zero,
                         const int max_subdivision)
{
  solver_info->max_itr_stag      = max_itr_stag;
  solver_info->max_itr_hardening = max_itr_hardening;
  solver_info->max_itr_M         = max_itr_M;
  solver_info->tol_hardening     = tol_hardening;
  solver_info->tol_M             = tol_M;
  solver_info->computer_zero     = computer_zero;
  solver_info->max_subdivision   = max_subdivision;
  solver_info->debug             = false;
}

/// print GCM solver info (numerical parameters)
///
/// \param[in] solver_info GcmSolverInfo structure
void print_gcm_solver_info(const GcmSolverInfo *solver_info)
{
  printf("-----------------------------------------------------------\n");
  printf("crystal plasticity solver info\n");
  printf("-----------------------------------------------------------\n");
  printf("max_itr_stag      = %d\n", solver_info->max_itr_stag);
  printf("max_itr_hardening = %d\n", solver_info->max_itr_hardening);
  printf("max_itr_M         = %d\n", solver_info->max_itr_M);
  printf("tol_hardening     = %e\n", solver_info->tol_hardening);
  printf("tol_M             = %e\n", solver_info->tol_M);
  printf("computer_zero     = %e\n", solver_info->computer_zero);
  printf("max_subdivision   = %d\n", solver_info->max_subdivision);
}