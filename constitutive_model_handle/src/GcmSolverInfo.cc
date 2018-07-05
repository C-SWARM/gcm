#include"GcmSolverInfo.h"

void set_gcm_solver_info(GcmSolverInfo *solver_info,
                        const int max_itr_stag, 
                        const int max_itr_hardening, 
                        const int max_itr_M,
                        const double tol_hardening, 
                        const double tol_M,
                        const double computer_zero)
{
  solver_info->max_itr_stag      = max_itr_stag;
  solver_info->max_itr_hardening = max_itr_hardening;
  solver_info->max_itr_M         = max_itr_M;
  solver_info->tol_hardening     = tol_hardening;
  solver_info->tol_M             = tol_M;
  solver_info->computer_zero     = computer_zero;
  solver_info->max_subdivision   = -1;
  solver_info->debug             = false;
}
