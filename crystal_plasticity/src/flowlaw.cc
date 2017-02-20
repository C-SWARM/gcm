#include "constitutive_model.h"
#include "flowlaw.h"
#include "material_properties.h"
#include <ttl/ttl.h>

#include <math.h>

namespace {
  template<int R, int D = 3, class S = double>
  using Tensor = ttl::Tensor<R, D, S>;

  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
}

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

  Tensor<2, 3, double*> C(C_in), S(S_in);
  Tensor<2> CS;   
  CS(i,j) = C(i,k) * S(k,j);   
  for(int a=0; a<slip->N_SYS; a++)
  {
    Tensor<2, 3, double*> P((slip->p_sys) + a*DIM_3x3);
    taus[a] = CS(i,j) * P(i,j);   //taus[a] equals the dot product of CS and P
  }
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
