#ifndef H__H__GCM_SOLVER_INFO__H__H
#define H__H__GCM_SOLVER_INFO__H__H

typedef struct
{
  int max_itr_stag;      // maximum number of iteration for staggerd NR
  int max_itr_hardening; // maximum number of iteration for hardening NR
  int max_itr_M;         // maximum number of iteration for M NR where (M = pFn_I)
  double tol_hardening;
  double tol_M;
  double computer_zero;
  double dt;
  int max_subdivision;
  bool debug;
} GcmSolverInfo;

void set_gcm_solver_info(GcmSolverInfo *solver_info,
                        const int max_itr_stag, 
                        const int max_itr_hardening, 
                        const int max_itr_M,
                        const double tol_hardening, 
                        const double tol_M, 
                        const double computer_zero,
                        const double dt);
#endif
