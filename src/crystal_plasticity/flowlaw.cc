#include "constitutive_model.h"
#include "flowlaw.h"
#include "material_properties.h"

#include <math.h>


// input double *gamma_dots, int N_SYS
// output double *gamma_dot
double compute_gamma_dot(double *gamma_dots, int N_SYS)
{
  int err = 0; 
  double gamma_dot = 0.0;
  for(int a = 0; a<N_SYS; a++)
    gamma_dot += fabs(gamma_dots[a]);

  return gamma_dot;  
}

int compute_gamma_dots(double *gamma_dots, double *taus, double g, 
                       MATERIAL_CRYSTAL_PLASTICITY *mat_p)
{
  int err = 0;
  double mm = 1.0/mat_p->m - 1.0;
  for(int a = 0; a<mat_p->slip->N_SYS; a++)
  {
    double A = taus[a]/g;
    gamma_dots[a] = mat_p->gamma_dot_0*A*pow(fabs(A), mm);
  }
  return err;
}       

int compute_tau_alphas(double *taus, double *C_in, double *S_in, SLIP_SYSTEM *slip)
{
  int err = 0;

  Matrix(double) CS,C,S,P;
  Matrix_construct_redim(double,CS,DIM_3,DIM_3);
  C.m_row = C.m_col = DIM_3; C.m_pdata = C_in;
  S.m_row = S.m_col = DIM_3; S.m_pdata = S_in;  
  P.m_row = P.m_col = DIM_3;   
  Matrix_AxB(CS,1.0,0.0,C,0,S,0);    
  for(int a=0; a<slip->N_SYS; a++)
  {
    P.m_pdata = (slip->p_sys) + a*DIM_3x3;
    Matrix_ddot(CS,P,taus[a]);
  }
  Matrix_cleanup(CS);
  return err;
}

int compute_d_gamma_d_tau(double *dgamma_dtaus, double g, double *taus, MATERIAL_CRYSTAL_PLASTICITY *mat_p)
{
  int err = 0;
  double mm = 1.0/mat_p->m - 1.0;
  for(int a = 0; a<mat_p->slip->N_SYS; a++)
    dgamma_dtaus[a] = mat_p->gamma_dot_0/mat_p->m/g*pow(fabs(taus[a]/g), mm);

  return err;
}
