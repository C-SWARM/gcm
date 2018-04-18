/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  Aaron Howell, [1], <ahowell3@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN

#include "constitutive_model.h"
#include "crystal_plasticity_integration.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "flowlaw.h"
#include "hardening.h"
#include "construct_linearization_parts.h"

#define D_GAMMA_D 0.005
#define D_GAMMA_TOL 1.25

#include <math.h>

long ODE_ops_EXA_metric = 0;

int set_crystal_plasticity_solver_info(CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info,
                                       int max_itr_stag, int max_itr_hardening, int max_itr_M,
                                       double tol_hardening, double tol_M, double computer_zero)
{
  int err = 0;
  solver_info->max_itr_stag      = max_itr_stag;
  solver_info->max_itr_hardening = max_itr_hardening;
  solver_info->max_itr_M         = max_itr_M;
  solver_info->tol_hardening     = tol_hardening;
  solver_info->tol_M             = tol_M;
  solver_info->computer_zero     = computer_zero;
  solver_info->max_subdivision   = -1;
  solver_info->debug             = 0;
  return err;
}

int print_crystal_plasticity_solver_info(CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info)
{
  int err = 0;
  printf("-----------------------------------------------------------\n");
  printf("crystal plasticity solver info\n");
  printf("-----------------------------------------------------------\n");
  printf("max_itr_stag      = %d\n", solver_info->max_itr_stag);
  printf("max_itr_hardening = %d\n", solver_info->max_itr_hardening);
  printf("max_itr_M         = %d\n", solver_info->max_itr_M);
  printf("tol_hardening     = %e\n", solver_info->tol_hardening);
  printf("tol_M             = %e\n", solver_info->tol_M);
  printf("computer_zero     = %e\n", solver_info->computer_zero);
  return err;
}

/// \param[out] R Vector with 10 elements
/// \param[in] M_in
/// \param[in] MI_in
/// \param[in] pFn_in
/// \param[in] pFnI_in
/// \param[in] N_in
/// \param[in] A_in
/// \param[in] gamma_dots
/// \param[in] lambda
/// \param[in] dt
/// \param[in] slip
/// \return non-zero on internal error
int compute_residual_M(double *R_out, double *M_in, double *MI_in, double *pFn_in, double *pFnI_in,
                       double *N_in, double *A_in,
                       double *gamma_dots,double lambda,double dt,
                       SLIP_SYSTEM *slip)
{
  int err = 0;

  Tensor<1, 10, double*> R(R_out);
  TensorA<2> M(M_in);
  TensorA<2> MI(MI_in);
  TensorA<2> pFn(pFn_in);
  TensorA<2> pFnI(pFnI_in);
  TensorA<2> N(N_in);
  TensorA<2> A(A_in);

  Tensor<2> Ru={},MIT,MIpFnN,Ident,AM,ATRu;

  for(int a=0; a<slip->N_SYS; a++)
  {
    TensorA<2> P((slip->p_sys) + a*DIM_3x3);
    Ru(i,j) += gamma_dots[a]*P(i,j);
  }

  MIT(i,j) = MI(j,i).to(i,j);
  MIpFnN(i,j) = MI(i,k) * pFn(k,l) * N(l,j);
  double det_MIpFnN = ttl::det(MIpFnN);
  Ident(i,j) = ttl::identity(i,j);
  AM(i,j) = A(i,k)*M(k,j);
  Ru(i,j) = dt*Ru(i,j) + AM(i,j) - Ident(i,j);
  ATRu(i,j) = A(k,i).to(i,k)*Ru(k,j);

  for(int a=0; a<DIM_3x3; a++)
    R.data[a] = ATRu.get(a) - lambda*det_MIpFnN*MIT.get(a);

  R.data[DIM_3x3] = det_MIpFnN - 1.0;

  return err;
}

/// \param[out] K_out
/// \param[in] M_in
/// \param[in] MI_in
/// \param[in] pFn_in
/// \param[in] eFnp1_in
/// \param[in] Fa_in
/// \param[in] C_in
/// \param[in] N_in
/// \param[in] A_in
/// \param[in] drdtaus
/// \param[in] lambda
/// \param[in] dt
/// \param[in] elasticity
/// \param[in] slip
/// \return non-zero on internal error
int construct_tangent_M(double *K_out, double *M_in, double *MI_in, double *pFn_in, double *eFnp1_in, double *Fa_in, double *C_in,
                        double *N_in, double *A_in, double *drdtaus,double lambda, double dt,
                        ELASTICITY *elasticity,SLIP_SYSTEM *slip)
{
  int err = 0;

  TensorA<2> MI(MI_in), pFn(pFn_in), eFnp1(eFnp1_in), Fa(Fa_in), N(N_in), A(A_in);

  Tensor<4> Kuu = {};
  Tensor<4> Kuu_a = {};
  Tensor<4> Kuu_b = {};

  Tensor<2> AA,det_MIpFnN_MI,MIpFnN;
  AA(i,j) = eFnp1(k,i).to(i,k)*Fa(k,j);

  // compute dt*sum dgamm/dtau Dtau[dM]*Pa
  for(int a=0; a<slip->N_SYS; a++)
  {
    compute_Kuu_a(Kuu_a.data,AA.data,elasticity->S,C_in,(slip->p_sys)+a*DIM_3x3,elasticity->L,drdtaus[a]);
    Kuu(i,j,k,l) += dt * Kuu_a(i,j,k,l);
  }

  // compute AT*(dt*sum dgamm/dtau Dtau[dM]*Pa + A*dM
  Kuu_a(i,j,k,l) = A(i,m).to(m,i) * (Kuu(m,j,k,l) + A(m,k) * ttl::identity(j,l));

  MIpFnN(i,j) = MI(i,k) * pFn(k,l) * N(l,j);
  double det_MIpFnN = ttl::det(MIpFnN);

  // compute tr(MI*dM)*MI' + (MI*dU*MI)'
  err += compute_Kuu_b(Kuu_b.data,MI.data);

  Tensor<2, 10, double*> K(K_out);

  for(int a=0; a<DIM_3x3; a++)
  {
    for(int b=0; b<DIM_3x3; b++)
      K(a,b) = Kuu_a.get(a*DIM_3x3+b) + lambda * det_MIpFnN * Kuu_b.get(a*DIM_3x3+b);
  }

  det_MIpFnN_MI(i,j) = det_MIpFnN * MI(j,i).to(i,j);
  for(int a=0; a<DIM_3x3; a++)
  {
    K(a,DIM_3x3) = -det_MIpFnN_MI.get(a);  //sets last row and column to -det_MIpFn_MI
    K(DIM_3x3,a) = -det_MIpFnN_MI.get(a);
  }

  K(DIM_3x3,DIM_3x3) = 0.0;   //K(9,9) = 0.0

  return err;
}

/// \param[out] M_out
/// \param[in] lambda
/// \param[in] pFn_in
/// \param[in] pFnI_in
/// \param[in] Fa_in
/// \param[in] N_in
/// \param[in] A_in
/// \param[in] g_np1_k
/// \param[in] dt
/// \param[in] mat
/// \param[in] elasticity
/// \param[in] solver_info
/// \param[in] is_it_cnvg
/// \param[in] norm_R_0
/// \param[in] d_gamma
/// \param[in] itr_stag
/// \return non-zero on internal error
double Newton_Rapson4M(double *M_out, double *lambda,
                    double *pFn_in, double *pFnI_in, double *Fa_in,
                    double *N_in, double *A_in,
                    double g_np1_k, double dt,
                    MATERIAL_CONSTITUTIVE_MODEL *mat,
                    ELASTICITY *elasticity,
                    CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info,
                    int *is_it_cnvg,
                    double *norm_R_0,
                    double *d_gamma,
                    int itr_stag)
{
  int err = 0;
  *is_it_cnvg = 0;

  SLIP_SYSTEM *slip = mat->mat_p->slip;

  TensorA<2> M(M_out);
  TensorA<2> Fa(Fa_in);

  Tensor<2> eFnp1,MI,C;

  double *taus         = new double[slip->N_SYS];
  double *gamma_dots   = new double[slip->N_SYS];
  double *dgamma_dtaus = new double[slip->N_SYS];

  Tensor<2, 10, double> K = {};
  Tensor<2, 10, double> KI = {};
  Tensor<1, 10, double> R = {};
  Tensor<1, 10, double> du = {};

  double norm_R = 0.0;
  double eng_norm_0 = 0.0; // energy norm
  // energy norm is used to check convergence
  // eng_norm = fabs(R*du)
  //
  // This check requires when variation of pF is very small
  // because when initial error (norm_R_0) is very small,
  // relative error (norm_R/norm_R_0) becomes large and Newton Rapson
  // is not convered.

  int is_it_cnvg_on_eng_norm = 0;

  int cnt = 0;
  double norm_R_n;

  for(int a = 0; a<solver_info->max_itr_M; a++)
  {
    cnt++;
    eFnp1(i,j) = Fa(i,k) * M(k,j);

    err += inv(M, MI);

    C(i,j) = eFnp1(k,i).to(i,k) * eFnp1(k,j);
    err += elasticity->update_elasticity(elasticity, eFnp1.data, 1); // compute stiffness also
    if(err != 0 )
      break;

    Tensor<2, 3, double*> S(elasticity->S);

    // compute taus
    err += compute_tau_alphas(taus,C.data, S.data, slip);
    // taus -> gamma_dots
    err += compute_gamma_dots(gamma_dots, taus, g_np1_k, mat->mat_p);
    // gamma_dots,taus -> dgamma_dtaus
    err += compute_d_gamma_d_tau(dgamma_dtaus, g_np1_k, taus, mat->mat_p);

    // compute R (risidual)
    err += compute_residual_M(R.data,M.data,MI.data,pFn_in,pFnI_in,
                              N_in, A_in, gamma_dots, *lambda, dt, slip);
    *d_gamma = 0.0;

    for(int s = 0; s<slip->N_SYS; s++)
      (*d_gamma) += dt*fabs(gamma_dots[s]);

    norm_R = R(i) * R(i);   // norm_R equals the dot product of R and R
    norm_R = sqrt(norm_R);

    if(a==0)
      norm_R_n = norm_R;

    if(a==0 && itr_stag==0)
    {
      *norm_R_0 = norm_R;

      if(*norm_R_0<solver_info->computer_zero)
        *norm_R_0 = solver_info->computer_zero;
    }

    if(norm_R/(*norm_R_0) < solver_info->tol_M)
      break;

    if(norm_R > norm_R_n)
    {
      if(DEBUG_PRINT_STAT)
        printf("Integration algorithm is diverging (M: %e -> %e)\n", norm_R_n, norm_R);
      break;
    }

    if((*d_gamma)/D_GAMMA_D>D_GAMMA_TOL && solver_info->max_subdivision>1)
      break;

    // compute dR/dM
    err += construct_tangent_M(K.data, M_out, MI.data,pFn_in,eFnp1.data,Fa_in,C.data,
                               N_in, A_in, dgamma_dtaus,*lambda, dt,elasticity,slip);

    // solve system of equation

    try
    {
      KI(i,j) = ttl::inverse(K)(i,j);  // attempt to take the inverse
    }
    catch(const int inverseException)  // no inverse exists
    {
      if(DEBUG_PRINT_STAT)
        printf( "Matrix is singular. The solution (Crystal plasticity) could not be computed.\n");

      err += 1;
      break;
    }

    du(i) = -KI(i,j) * R(j);

    // update M and compute energy norm = fabs(R*du)
    double eng_norm = 0.0;
    for(int b=0; b<DIM_3x3; b++)
    {
      M.data[b] += du.data[b];
      eng_norm += du.data[b]*R.data[b];
    }
    *lambda += du.data[DIM_3x3];
    eng_norm += du.data[DIM_3x3]*R.data[DIM_3x3];
    eng_norm = fabs(eng_norm);

    // set initial energy norm
    if(a==0)
      eng_norm_0 = eng_norm;

    if(eng_norm_0<solver_info->computer_zero)
      eng_norm_0 = solver_info->computer_zero;

    // check convergence based on energy norm
    if(eng_norm/eng_norm_0 < (solver_info->tol_M)*(solver_info->tol_M))
    {
      if(DEBUG_PRINT_STAT)
        printf("converge on energe norm %e\n", eng_norm/eng_norm_0);

      is_it_cnvg_on_eng_norm = 1;
      break;

    }
    norm_R_n = norm_R;
  }

  if(DEBUG_PRINT_STAT)
  {
    printf("itr = %d\n", cnt-1);
    printf("residual: |R| = %e, |R0| = %e, |R/R0| = %e\n", norm_R, *norm_R_0, norm_R/(*norm_R_0));
  }

  if(((norm_R/(*norm_R_0) < solver_info->tol_M) && (err==0)) || is_it_cnvg_on_eng_norm)
    *is_it_cnvg = 1;

  delete [] taus;
  delete [] gamma_dots;
  delete [] dgamma_dtaus;

  return err;
}

int Newton_Rapson_hardening(double *g_np1, double *eFnp1_in,
                            double g_n, double g_np1_k, double dt,
                            MATERIAL_CONSTITUTIVE_MODEL *mat,
                            ELASTICITY *elasticity,
                            CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info)
{
  int err = 0;

  SLIP_SYSTEM *slip = mat->mat_p->slip;

  // this will use pointers, no need Matrix_cleanup -->
  TensorA<2> eFnp1(eFnp1_in);
  // <-- this will use pointers, no need Matrix_cleanup

  Tensor<2> C;
  C = eFnp1(k, i)*eFnp1(k,j);

  double *taus         = new double[slip->N_SYS];
  double *gamma_dots   = new double[slip->N_SYS];


  // compute stress -->
  elasticity->update_elasticity(elasticity,eFnp1.data, 0);
  // <-- compute stress

  // compute taus
  err += compute_tau_alphas(taus,C.data, elasticity->S, slip);
  // taus -> gamma_dots
  err += compute_gamma_dots(gamma_dots, taus, g_np1_k, mat->mat_p);
  // gamma_dots -> gamma_dot_np1
  double gamma_dot_np1 =  compute_gamma_dot(gamma_dots, slip->N_SYS);
  // compute
  double gs_np1 = compute_gs_np1(gamma_dot_np1, mat->mat_p);
  // compute g_np1
  *g_np1  = compute_g_np1(gs_np1, g_n, dt, gamma_dot_np1, mat->mat_p);

  delete [] taus;
  delete [] gamma_dots;

  return err;
}

int staggered_Newton_Rapson_compute(double *pFnp1_out, double *M_out, double *g_out, double *lambda,
                            double *pFn_in, double *Fn_in, double *Fnp1_in, double *hFn_in, double *hFnp1_in,
                            double g_n, double dt,
                            MATERIAL_CONSTITUTIVE_MODEL *mat,
                            ELASTICITY *elasticity,
                            CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info,
                            double *d_gamma,
                            int *is_it_cnvg)
{
  int err = 0;

  ++ODE_ops_EXA_metric;
  
  double g_np1_k   = g_n;
  double g_np1_kp1 = g_n;

  double g_0, err_g = 0.0, err_g_n = 0.0;
  g_0 = (mat->mat_p)->g0;

  TensorA<2> pFnp1(pFnp1_out), M(M_out), pFn(pFn_in), Fn(Fn_in), Fnp1(Fnp1_in), hFn(hFn_in), hFnp1(hFnp1_in);
  Tensor<2> eFn,eFnI,pFnI,FnI,Fr,FrI,Fa,MI,hFnI,eFnp1,hFnp1I,A,N,NI;

  err += inv(pFn,pFnI);
  err += inv(hFn,hFnI);
  err += inv(Fn, FnI);

  // compute eFn
  eFn = Fn(i,k)*hFnI(k,l)*pFnI(l,j);
  Fr = Fnp1(i,k)*FnI(k,j);

  err += inv(Fr, FrI);
  // guess initial M
  err += inv(eFn, eFnI);
  M = eFnI(i,k)*FrI(k,l)*eFn(l,j);
  // compute Fa
  Fa = Fr(i,k)*eFn(k,j);
  // compute N : N = hFn*hFnp1I
  err += inv(hFnp1, hFnp1I);
  N = hFn(i,k)*hFnp1I(k,j);
  // compute A : A = pFn*N_I*pFn_I
  err += inv(N, NI);
  A = pFn(i,k)*NI(k,l)*pFnI(l,j);

  int is_it_cnvg_M = 0;
  double norm_R_0 = 0.0;

  // if there is err, do nothing
  int max_itr = solver_info->max_itr_stag;
  if(err>0)
    max_itr = -1;

  for(int k=0; k<max_itr; k++)
  {
    if(DEBUG_PRINT_STAT)
      printf("staggered iter = %d:\n", k);

    if(k==0)
    {
      // initial guess M is computed based on assumption that eFnp1 = eFn
      // so that initial hardening is computed using eFn
      err += Newton_Rapson_hardening(&g_np1_kp1, eFn.data,
                            g_n, g_np1_k, dt, mat, elasticity, solver_info);

      err_g = sqrt((g_np1_kp1 - g_np1_k)*(g_np1_kp1 - g_np1_k));
      g_np1_k = g_np1_kp1;
    }
    err += Newton_Rapson4M(M.data, lambda,
                           pFn.data, pFnI.data, Fa.data,N.data,A.data,
                           g_np1_k, dt, mat, elasticity, solver_info,&is_it_cnvg_M,
                           &norm_R_0,d_gamma,k);

    if(err>0)
      break;

    // compute eFnp1 = Fa*M
    eFnp1 = Fa(i,r)*M(r,j);
    err += Newton_Rapson_hardening(&g_np1_kp1, eFnp1.data,
                            g_n, g_np1_k, dt, mat, elasticity, solver_info);

    err_g = sqrt((g_np1_kp1 - g_np1_k)*(g_np1_kp1 - g_np1_k));
    g_np1_k = g_np1_kp1;

    if(k==0 && err_g>solver_info->computer_zero)
      err_g_n = err_g;

    if(DEBUG_PRINT_STAT)
    {
      printf("err(g_np1) \t= %e\n", err_g);
      printf("w_max \t= %e\n", (*d_gamma)/D_GAMMA_D);
      printf("\n");
    }

    if((err_g/g_0)<solver_info->tol_hardening && is_it_cnvg_M)
      break;

    if(err_g > err_g_n)
    {
      if(DEBUG_PRINT_STAT)
        printf("Integration algorithm is diverging (g: %e -> %e)\n", err_g_n, err_g);
      break;
    }
    if((*d_gamma)/D_GAMMA_D > D_GAMMA_TOL && solver_info->max_subdivision > 1)
    {
      if(DEBUG_PRINT_STAT)
        printf("Need smaller dt w_max = %e\n", (*d_gamma)/D_GAMMA_D);
      break;
    }

    err_g_n = err_g;
  }

  if((err_g/g_0)<solver_info->tol_hardening && is_it_cnvg_M)
  {
    err += inv(M,MI);
    pFnp1 = MI(i,k)*pFn(k,l)*N(l,j);
    *g_out = g_np1_kp1;
    *is_it_cnvg = 1;
  }
  else
    err += 1;

  return err;
}

int staggered_Newton_Rapson_subdivision(double *pFnp1_out, double *M_out, double *g_out, double *lambda,
                            double *pFn_in, double *Fn_in, double *Fnp1_in, double *hFn_in, double *hFnp1_in,
                            double g_n, double dt,
                            MATERIAL_CONSTITUTIVE_MODEL *mat,
                            ELASTICITY *elasticity,
                            CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info,
                            int stepno,
                            double *d_gamma,
                            int *is_it_cnvg)
{
  int err = 0;

  Tensor<2> dF,F_s,Fn_s,pFn_s;

  double dt_s = dt/stepno;

  for(int b=0; b<DIM_3x3; b++)
  {
       dF.data[b] = (Fnp1_in[b] - Fn_in[b])/stepno;
     Fn_s.data[b] =  Fn_in[b];
    pFn_s.data[b] = pFn_in[b];
  }

  double g_n_s = g_n;

  int is_sub_cnvg = 0;

  for(int a=0; a<stepno; a++)
  {
    if(DEBUG_PRINT_STAT)
      printf("-----------------------\n sub step: (%d/%d) \n-----------------------\n", a+1, stepno);

    is_sub_cnvg = 0;
    F_s = Fn_s(i,j) + dF(i,j);

    err += staggered_Newton_Rapson_compute(pFnp1_out, M_out, g_out,lambda,
                                         pFn_s.data,
                                          Fn_s.data,
                                           F_s.data,
                                         hFn_in, hFnp1_in,
                                         g_n_s,dt_s,mat,
                                         elasticity,solver_info,d_gamma,&is_sub_cnvg);
    if((*d_gamma)/D_GAMMA_D>D_GAMMA_TOL || !is_sub_cnvg)
    {
      err = 1;
      break;
    }

    for(int b=0; b<DIM_3x3; b++)
    {
      pFn_s.data[b] = pFnp1_out[b];
       Fn_s.data[b] = F_s.data[b];
    }
    g_n_s = *g_out;
  }

  *is_it_cnvg = is_sub_cnvg;

  return err;
}

/// integrate crystal plasticity from pFn to pFnp1 for general cases (includes thermal expansions)
///
/// In integrating plastic part of deformation gradient and hardening,
/// Both non-linear iterations are staggered until two rediduals (M = pFr_I and g = hardening) are converged.
/// If one of the iterations is diverging, subdivision steps can be applied if solver_info include
/// subdivision number greater than 1. In the subdivision steps, loading (Fnp1) is linearly factorized as smaller
/// time step size used. The subdivision step size increases 2, 4, 8, 18, ... upto solver_infor.max_subdivision.
///
/// \param[out] pFnp1_out deformation gradient (plastic) at time step = n + 1
/// \param[out] M_out deformation gradient (pFr_I) at time step = n + 1
/// \param[out] g_out updated hardening at time step =  n + 1
/// \param[out] lambda updated Lagrange multiplier at time step =  n + 1
/// \param[in] pFn_in deformation gradient (plastic) at time step = n
/// \param[in] Fn_in deformation gradient (total) at time step = n
/// \param[in] Fnp1_in deformation gradient (total) at time step = n + 1
/// \param[in] hFn_in deformation gradient (due to thermal expansion) at time step = n
/// \param[in] hFnp1_in deformation gradient (due to thermal expansion) at time step = n + 1
/// \param[in] g_n hardening at time step =  n
/// \param[in] dt time step size
/// \param[in] mat material parameters
/// \param[in] elasticity object to compute elasticity (stress)
/// \param[in] solver_info defines numerical parameters
/// \return non-zero on interal error
int staggered_Newton_Rapson_generalized(double *pFnp1_out, double *M_out, double *g_out, double *lambda,
                                        double *pFn_in, double *Fn_in, double *Fnp1_in, double *hFn_in, double *hFnp1_in,
                                        double g_n, double dt,
                                        MATERIAL_CONSTITUTIVE_MODEL *mat,
                                        ELASTICITY *elasticity,
                                        CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info)
{
  int err = 0;
  double lambda_in = *lambda;
  int is_it_cnvg = 0;
  double d_gamma;
  err += staggered_Newton_Rapson_compute(pFnp1_out, M_out, g_out,lambda,
                                         pFn_in,Fn_in,Fnp1_in,hFn_in,hFnp1_in,g_n,dt,mat,
                                         elasticity,solver_info,&d_gamma,&is_it_cnvg);

  if(solver_info->max_subdivision<2)
  {
    if(is_it_cnvg==0)
      printf("Integration algorithm is diverging\n");
    return err;
  }

  double w = (d_gamma)/D_GAMMA_D;
  if(w < D_GAMMA_TOL && is_it_cnvg)
    return err;

  if(!is_it_cnvg || w>D_GAMMA_TOL)
  {
    if(DEBUG_PRINT_STAT)
    {
      printf("WARNING: Time stepping size is adjusted in Crystal plasticity integration (w=%e)\n", w);
      printf("Enterring sub division ..\n");
    }

    int stepno = 2;
    while(w>D_GAMMA_TOL || !is_it_cnvg)
    {
      *lambda = lambda_in;
      err = 0;
      is_it_cnvg = 0;
      err += staggered_Newton_Rapson_subdivision(pFnp1_out, M_out, g_out,lambda,
                                                 pFn_in,Fn_in,Fnp1_in,hFn_in,hFnp1_in,g_n,dt,
                                                 mat,elasticity,solver_info,stepno,&d_gamma,&is_it_cnvg);
      stepno *= 2;
      w = d_gamma/D_GAMMA_D;

      if(stepno>(solver_info->max_subdivision))
      {
        if(w>D_GAMMA_TOL)
        {
          if(DEBUG_PRINT_STAT)
          {
            printf("WARNING: reach maximum number of sub division\n");
            printf("Integration algorithm updates values computed at the final step\n");
          }
          err++;
        }
        break;
      }
    }
  }

  if(is_it_cnvg==0)
    printf("After subdivision: Integration algorithm is diverging\n");

  return err;
}

int Fnp1_Implicit(double *Fnp1_out, double *Fn_in, double *L_in, double dt)
{
  int err = 0;
  double I[9];
  I[0] = I[4] = I[8] = 1.0;
  I[1] = I[2] = I[3] = I[5] = I[6] = I[7] = 0.0;

  for(int a=0; a<DIM_3x3; a++)
    I[a] -= dt*L_in[a];

  TensorA<2> A(I), Fnp1(Fnp1_out), Fn(Fn_in);
  Tensor<2> AI;

  err += inv(A, AI);
  Fnp1 = AI(i,k)*Fn(k,j);
  return err;
}

/// integrate crystal plasticity from pFn to pFnp1
///
/// Integrating plastic part of deformation gradient and hardening
/// same as staggered_Newton_Rapson_generalized except not includs thermal expansions.
///
/// \param[out] pFnp1_out deformation gradient (plastic) at time step = n + 1
/// \param[out] M_out deformation gradient (pFr_I) at time step = n + 1
/// \param[out] g_out updated hardening at time step =  n + 1
/// \param[out] lambda updated Lagrange multiplier at time step =  n + 1
/// \param[in] pFn_in deformation gradient (plastic) at time step = n
/// \param[in] Fn_in deformation gradient (total) at time step = n
/// \param[in] Fnp1_in deformation gradient (total) at time step = n + 1
/// \param[in] g_n hardening at time step =  n
/// \param[in] dt time step size
/// \param[in] mat material parameters
/// \param[in] elasticity object to compute elasticity (stress)
/// \param[in] solver_info defines numerical parameters
/// \return non-zero on interal error
int staggered_Newton_Rapson(double *pFnp1_out, double *M_out, double *g_out, double *lambda,
                            double *pFn_in, double *Fn_in, double *Fnp1_in,
                            double g_n, double dt,
                            MATERIAL_CONSTITUTIVE_MODEL *mat,
                            ELASTICITY *elasticity,
                            CRYSTAL_PLASTICITY_SOLVER_INFO *solver_info)
{
  double hFn[9]; // hFn*hFnp1_I
  hFn[1] = hFn[2] = hFn[3] = hFn[5] = hFn[6] = hFn[7] = 0.0;
  hFn[0] = hFn[4] = hFn[8] = 1.0;

  double hFnp1[9]; // hFn*hFnp1_I
  hFnp1[1] = hFnp1[2] = hFnp1[3] = hFnp1[5] = hFnp1[6] = hFnp1[7] = 0.0;
  hFnp1[0] = hFnp1[4] = hFnp1[8] = 1.0;

  int err = staggered_Newton_Rapson_generalized(pFnp1_out, M_out, g_out,lambda,
                                                pFn_in, Fn_in, Fnp1_in, hFn, hFnp1,
                                                g_n, dt, mat, elasticity, solver_info);
  return err;
}
