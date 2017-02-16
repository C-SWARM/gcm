#include "constitutive_model.h"
#include "crystal_plasticity_integration.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "flowlaw.h"
#include "hardening.h"
#include "construct_linearization_parts.h"
#include <ttl/ttl.h>

#define D_GAMMA_D 0.005
#define D_GAMMA_TOL 1.25

#include <math.h>

namespace {
  template<int R, int D = 3, class S = double>
  using Tensor = ttl::Tensor<R, D, S>;

  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;
}

/// compute C = A*B for 3 by 3 matrix
///
/// \param[out] C [3 by 3] matrix, result of A*B
/// \param[in]  A [3 by 3] input matrix to be multiplied
/// \param[in]  B [3 by 3] input matrix to be multiplied
/// \return non-zero on interal error
int Matrix_AxB_3by3(Matrix(double) *C, const Matrix(double) *A, Matrix(double) *B)
{
  int err = 0;
  for(int ia=1; ia<=DIM_3; ia++)
  {
    for(int ib=1; ib<=DIM_3; ib++)
    {
      Mat_v(*C, ia, ib) = 0.0;
      for(int ic=1; ic<=DIM_3; ic++)
        Mat_v(*C, ia, ib) += Mat_v(*A, ia, ic)*Mat_v(*B, ic, ib);
    }
  }
  return err;
}

/// compute C = A'*B for 3 by 3 matrix
///
/// \param[out] C [3 by 3] matrix, result of A*B
/// \param[in]  A [3 by 3] input matrix to be multiplied
/// \param[in]  B [3 by 3] input matrix to be multiplied
/// \return non-zero on interal error
int Matrix_ATxB_3by3(Matrix(double) *C, const Matrix(double) *A, Matrix(double) *B)
{
  int err = 0;
  for(int ia=1; ia<=DIM_3; ia++)
  {
    for(int ib=1; ib<=DIM_3; ib++)
    {
      Mat_v(*C, ia, ib) = 0.0;
      for(int ic=1; ic<=DIM_3; ic++)
        Mat_v(*C, ia, ib) += Mat_v(*A, ic, ia)*Mat_v(*B, ic, ib);
    }
  }
  return err;
}

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
                        
int compute_residual_M(double *R, double *M_in, double *MI_in, double *pFn_in, double *pFnI_in, 
                       double *N_in, double *A_in,
                       double *gamma_dots,double lambda,double dt,
                       SLIP_SYSTEM *slip)
{
  int err = 0;
  
  Tensor<2, 3, double*> M(M_in);
  Tensor<2, 3, double*> MI(MI_in);
  Tensor<2, 3, double*> pFn(pFn_in);
  Tensor<2, 3, double*> pFnI(pFnI_in);
  Tensor<2, 3, double*> N(N_in);
  Tensor<2, 3, double*> A(A_in);     
    
  enum {Ru,MIT,MIpFnN,Ident,AM,ATRu,F2end};
  Tensor<2> F2[F2end];   //declare an array of tensors and initialize them to 0
  for (int a = 0; a < F2end; a++) {
    F2[a] = {};
  } 
  
  //Tensor<2> P = {};   //can try putting this inside the loop (Tensor<2,double*> P((slip.p_sys) + a*DIM_3x3);)
          
  for(int a=0; a<slip->N_SYS; a++)
  {
    Tensor<2, 3, double*> P((slip->p_sys) + a*DIM_3x3);  //declaring p was originally done outside the loop, may want to move it back out
    F2[Ru](i,j) += gamma_dots[a] * P(i,j);
  }

  F2[MIT](i,j) = MI(j,i).to(i,j);     
  F2[MIpFnN](i,j) = MI(i,k) * pFn(k,l) * N(l,j);
  //Matrix_Tns2_AxBxC(F2[MIpFnN],1.0,0.0,MI,pFn,N);
  double det_MIpFnN = ttl::det(F2[MIpFnN]); 
  F2[Ident](i,j) = ttl::identity(i,j);
  F2[AM](i,j) = A(i,k) * M(k,j);
  //Matrix_AxB_3by3(F2+AM,&A,&M); // <= Matrix_AxB(F2[AM], 1.0,0.0,A,0,M,0);
  
  F2[Ru](i,j) = dt*F2[Ru](i,j) + F2[AM](i,j) - F2[Ident](i,j);
/*
  for(int a=0; a<DIM_3x3; a++)
    F2[Ru].m_pdata[a] = dt*F2[Ru].m_pdata[a] + F2[AM].m_pdata[a] - F2[Ident].m_pdata[a];  
*/
  
    F2[ATRu](i,j) = A(k,i).to(i,k) * F2[Ru](k,j);
    //Matrix_ATxB_3by3(F2+ATRu,&A,F2+Ru); // <= Matrix_AxB(F2[ATRu], 1.0,0.0,A,1,F2[Ru],0);  
  
  for(int a=0; a<DIM_3x3; a++)
    R[a] = F2[ATRu].get(a) - lambda * det_MIpFnN * F2[MIT].get(a);
  

  R[DIM_3x3] = det_MIpFnN - 1.0;
          
  return err;
}

int construct_tangent_M(double *K_out, double *M_in, double *MI_in, double *pFn_in, double *eFnp1_in, double *Fa_in, double *C_in,
                        double *N_in, double *A_in, double *drdtaus,double lambda, double dt,
                        ELASTICITY *elasticity,SLIP_SYSTEM *slip)
{
  int err = 0;

  Matrix(double) M,MI,pFn,eFnp1,Fa,N,A;
      M.m_row =     M.m_col = DIM_3;     M.m_pdata =    M_in;
     MI.m_row =    MI.m_col = DIM_3;    MI.m_pdata =    MI_in;
    pFn.m_row =   pFn.m_col = DIM_3;   pFn.m_pdata =   pFn_in; 
  eFnp1.m_row = eFnp1.m_col = DIM_3; eFnp1.m_pdata = eFnp1_in;
     Fa.m_row =    Fa.m_col = DIM_3;    Fa.m_pdata =    Fa_in;     
      N.m_row =     N.m_col = DIM_3;     N.m_pdata = N_in;  // N = hFn*hFnp1_I
      A.m_row =     A.m_col = DIM_3;     A.m_pdata = A_in;  // A = pFn*N_I*pFn_I     
  Matrix(double) Kuu, Kuu_a, Kuu_b;
  Matrix_construct_init(double,Kuu,  DIM_3x3,DIM_3x3,0.0);
  Matrix_construct_init(double,Kuu_a,DIM_3x3,DIM_3x3,0.0); 
  Matrix_construct_init(double,Kuu_b,DIM_3x3,DIM_3x3,0.0);
  
  enum {AA,det_MIpFnN_MI,MIpFnN,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++)
    Matrix_construct_redim(double, F2[a],DIM_3,DIM_3);
  
  Matrix_ATxB_3by3(F2+AA,&eFnp1,&Fa); // <= Matrix_AxB(F2[AA],1.0,0.0,eFnp1,1,Fa,0);  
  
  // compute dt*sum dgamm/dtau Dtau[dM]*Pa 
  for(int a=0; a<slip->N_SYS; a++)
  {
    compute_Kuu_a(Kuu_a.m_pdata,F2[AA].m_pdata,elasticity->S,C_in,(slip->p_sys)+a*DIM_3x3,elasticity->L,drdtaus[a]); 
    Matrix_AplusB(Kuu,1.0,Kuu,dt,Kuu_a);
  }
  
  // compute AT*(dt*sum dgamm/dtau Dtau[dM]*Pa + A*dM  
  for(int ia=1; ia<=DIM_3; ia++)
  {
    for(int ic=1; ic<=DIM_3; ic++)
    {
      for(int id=1; id<=DIM_3; id++)
      {
        for(int ie=1; ie<=DIM_3; ie++)
        {
          Tns4_v(Kuu_a,ia,ic,id,ie) = 0.0;
          for(int ib=1; ib<=DIM_3; ib++)
            Tns4_v(Kuu_a,ia,ic,id,ie) += Mat_v(A,ib,ia)*(Tns4_v(Kuu,ib,ic,id,ie) + Mat_v(A,ib,id)*(ic==ie));
        }
      }
    }
  }  
  
  Matrix_Tns2_AxBxC(F2[MIpFnN],1.0,0.0,MI,pFn,N);
  double det_MIpFnN;
  Matrix_det(F2[MIpFnN], det_MIpFnN);  

  // compute tr(MI*dM)*MI' + (MI*dU*MI)'  
  err += compute_Kuu_b(Kuu_b.m_pdata,MI.m_pdata);  
  
  Matrix(double) K;
  K.m_row = K.m_col = DIM_3x3+1;  K.m_pdata = K_out;
  
  for(int a=1; a<=DIM_3x3; a++)
  {
    for(int b=1; b<=DIM_3x3; b++)
      Mat_v(K,a,b) = Mat_v(Kuu_a,a,b) + lambda*det_MIpFnN*Mat_v(Kuu_b,a,b);
  }
  
  Matrix_AeqBT(F2[det_MIpFnN_MI],det_MIpFnN,MI);
  for(int a=1; a<=DIM_3x3; a++)
  {
    Mat_v(K,a,DIM_3x3+1) = -F2[det_MIpFnN_MI].m_pdata[a-1];
    Mat_v(K,DIM_3x3+1,a) = -F2[det_MIpFnN_MI].m_pdata[a-1];    
  }

  Mat_v(K,DIM_3x3+1,DIM_3x3+1) = 0.0;  

  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);
  free(F2);  
  Matrix_cleanup(Kuu);
  Matrix_cleanup(Kuu_a);
  Matrix_cleanup(Kuu_b);  
  return err;
}

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
  
  // this will use pointers, no need Matrix_cleanup -->
  Matrix(double) M, Fa;
     M.m_row =    M.m_col = DIM_3;    M.m_pdata = M_out;   
    Fa.m_row =   Fa.m_col = DIM_3;   Fa.m_pdata =   Fa_in;
    
  // this will use pointers, no need Matrix_cleanup -->  
  
  enum {eFnp1,MI,C,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_redim(double, F2[a],DIM_3,DIM_3);
  } 
    
  enum {taus,gamma_dots,dgamma_dtaus,F1end};
  Matrix(double) *F1 = malloc(F1end*sizeof(Matrix(double)));
  for (int a = 0; a < F1end; a++) {
    Matrix_construct_redim(double, F1[a],slip->N_SYS, 1);
  } 
  
  Matrix(double) K, KI, R, du;
  Matrix_construct_redim(double,K, DIM_3x3+1,DIM_3x3+1);
  Matrix_construct_redim(double,KI,DIM_3x3+1,DIM_3x3+1);
  Matrix_construct_redim(double,R, DIM_3x3+1,1);
  Matrix_construct_redim(double,du,DIM_3x3+1,1);
  
  double norm_R;    
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
    Matrix_AxB_3by3(F2+eFnp1,&Fa,&M); // <= Matrix_AxB(F2[eFnp1],1.0,0.0,Fa,0,M,0);
          
    Matrix_inv(M, F2[MI]);
    Matrix_ATxB_3by3(F2+C,F2+eFnp1,F2+eFnp1); // <=  Matrix_AxB(F2[C],1.0,0.0,F2[eFnp1],1,F2[eFnp1],0);        
    elasticity->update_elasticity(elasticity,F2[eFnp1].m_pdata, 1); // compute stiffness also
    
    Matrix(double) S;
    S.m_col = S.m_row = 3; S.m_pdata = elasticity->S;

    // compute taus
    err += compute_tau_alphas(F1[taus].m_pdata,F2[C].m_pdata, S.m_pdata, slip);    
    // taus -> gamma_dots
    err += compute_gamma_dots(F1[gamma_dots].m_pdata, F1[taus].m_pdata, g_np1_k, mat->mat_p);
    // gamma_dots,taus -> dgamma_dtaus
    err += compute_d_gamma_d_tau(F1[dgamma_dtaus].m_pdata, g_np1_k, F1[taus].m_pdata, mat->mat_p);    

    // compute R (risidual)
    err += compute_residual_M(R.m_pdata,M.m_pdata,F2[MI].m_pdata,pFn_in,pFnI_in,
                              N_in, A_in, F1[gamma_dots].m_pdata, *lambda, dt, slip);
    *d_gamma = 0.0;
    
    for(int s = 0; s<slip->N_SYS; s++)
      (*d_gamma) += dt*fabs(F1[gamma_dots].m_pdata[s]);
                                          
    Matrix_ddot(R,R,norm_R);
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
    err += construct_tangent_M(K.m_pdata, M_out, F2[MI].m_pdata,pFn_in,F2[eFnp1].m_pdata,Fa_in,F2[C].m_pdata,
                               N_in, A_in, F1[dgamma_dtaus].m_pdata,*lambda, dt,elasticity,slip);
    
    // solve system of equation

    int info = 0;    

    Matrix_inv_check_err(K,KI,info);
    Matrix_AxB(du,-1.0,0.0,KI,0,R,0);
                    
    if(info <= 0)
    {
      // update M and compute energy norm = fabs(R*du)
      double eng_norm = 0.0;      
      for(int b=0; b<DIM_3x3; b++)
      {
        M.m_pdata[b] += du.m_pdata[b];
        eng_norm += du.m_pdata[b]*R.m_pdata[b];
      }
      *lambda += du.m_pdata[DIM_3x3];
      eng_norm += du.m_pdata[DIM_3x3]*R.m_pdata[DIM_3x3];
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
    else
    {
      if(DEBUG_PRINT_STAT)
        printf( "Matrix is singular. The solution (Crystal plasticity) could not be computed.\n");

      err += 1;
      break;
    }
  }
  
  if(DEBUG_PRINT_STAT)
  {  
    printf("itr = %d\n", cnt-1);
    printf("residual: |R| = %e, |R0| = %e, |R/R0| = %e\n", norm_R, *norm_R_0, norm_R/(*norm_R_0));    
  }
  
  if((norm_R/(*norm_R_0) < solver_info->tol_M) && err==0 || is_it_cnvg_on_eng_norm)
    *is_it_cnvg = 1;  
       
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);
  free(F2);   
     
  for(int a = 0; a < F1end; a++)
    Matrix_cleanup(F1[a]);
  free(F1);  
    
  Matrix_cleanup(K);
  Matrix_cleanup(KI);  
  Matrix_cleanup(R);
  Matrix_cleanup(du);
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
  Matrix(double) eFnp1; 
  eFnp1.m_row = eFnp1.m_col = DIM_3; eFnp1.m_pdata = eFnp1_in;
  // <-- this will use pointers, no need Matrix_cleanup
  
  Matrix(double) C;
  Matrix_construct_redim(double,C,DIM_3,DIM_3);

  enum {taus,gamma_dots,F1end};
  Matrix(double) *F1 = malloc(F1end*sizeof(Matrix(double)));
  for (int a = 0; a < F1end; a++) {
    Matrix_construct_redim(double, F1[a],slip->N_SYS, 1);
  } 
  
  // compute stress -->            
  Matrix_ATxB_3by3(&C,&eFnp1,&eFnp1); // <= Matrix_AxB(C, 1.0, 0.0, eFnp1,1,eFnp1,0);    
  elasticity->update_elasticity(elasticity,eFnp1.m_pdata, 0);
  // <-- compute stress
  
  // compute taus
  err += compute_tau_alphas(F1[taus].m_pdata,C.m_pdata, elasticity->S, slip);
  // taus -> gamma_dots
  err += compute_gamma_dots(F1[gamma_dots].m_pdata, F1[taus].m_pdata, g_np1_k, mat->mat_p);
  // gamma_dots -> gamma_dot_np1
  double gamma_dot_np1 =  compute_gamma_dot(F1[gamma_dots].m_pdata, slip->N_SYS);
  // compute 
  double gs_np1 = compute_gs_np1(gamma_dot_np1, mat->mat_p);
  // compute g_np1
  *g_np1  = compute_g_np1(gs_np1, g_n, dt, gamma_dot_np1, mat->mat_p);

  
  Matrix_cleanup(C);   
     
  for(int a = 0; a < F1end; a++)
    Matrix_cleanup(F1[a]);
  free(F1);
              
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
                  
  double g_np1_k   = g_n;
  double g_np1_kp1 = g_n;
  
  double g_0, err_g, err_g_n;
  g_0 = (mat->mat_p)->g0;

  // this will use pointers, no need Matrix_cleanup -->
  Matrix(double) pFnp1, M, pFn, Fn, Fnp1, hFn, hFnp1;
  pFnp1.m_row = pFnp1.m_col = DIM_3; pFnp1.m_pdata = pFnp1_out;   
      M.m_row =     M.m_col = DIM_3;     M.m_pdata = M_out;
    pFn.m_row =   pFn.m_col = DIM_3;   pFn.m_pdata = pFn_in;
     Fn.m_row =    Fn.m_col = DIM_3;    Fn.m_pdata = Fn_in;  
   Fnp1.m_row =  Fnp1.m_col = DIM_3;  Fnp1.m_pdata = Fnp1_in;
    hFn.m_row =   hFn.m_col = DIM_3;   hFn.m_pdata = hFn_in;
  hFnp1.m_row = hFnp1.m_col = DIM_3; hFnp1.m_pdata = hFnp1_in;        
  // <-- this will use pointers, no need Matrix_cleanup
  
  enum {eFn,eFnI,pFnI,FnI,Fr,FrI,Fa,MI,hFnI,eFnp1,hFnp1I,A,N,NI,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_redim(double, F2[a],DIM_3,DIM_3);
  } 
    
  // compute eFn
  int info = 0;
  Matrix_inv_check_err(pFn, F2[pFnI], info);
  err += info;
  Matrix_inv_check_err(hFn, F2[hFnI], info);
  err += info;  
  Matrix_Tns2_AxBxC(F2[eFn],1.0,0.0,Fn,F2[hFnI],F2[pFnI]);
  // compute Fr
  Matrix_inv_check_err(Fn, F2[FnI], info);
  err += info;
  Matrix_AxB_3by3(F2+Fr,&Fnp1,F2+FnI); // <= Matrix_AxB(F2[Fr],1.0,0.0,Fnp1,0,F2[FnI],0);
  Matrix_inv_check_err(F2[Fr], F2[FrI], info);
  err += info;
  // guess initial M
  Matrix_inv_check_err(F2[eFn], F2[eFnI], info);
  err += info;  
  Matrix_Tns2_AxBxC(M,1.0,0.0,F2[eFnI],F2[FrI],F2[eFn]);
  // compute Fa  
  Matrix_AxB_3by3(F2+Fa,F2+Fr,F2+eFn); // <= Matrix_AxB(F2[Fa],1.0,0.0,F2[Fr],0,F2[eFn],0);
  // compute N : N = hFn*hFnp1_I
  Matrix_inv_check_err(hFnp1, F2[hFnp1I], info);
  err += info;
  Matrix_AxB_3by3(F2+N,&hFn,F2+hFnp1I); // <= Matrix_AxB(F2[N],1.0,0.0,hFn,0,F2[hFnp1I],0);   
  // compute A : A = pFn*N_I*pFn_I   
  Matrix_inv_check_err(F2[N], F2[NI], info);
  err += info;
  Matrix_Tns2_AxBxC(F2[A],1.0,0.0,pFn,F2[NI],F2[pFnI]);    

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
      err += Newton_Rapson_hardening(&g_np1_kp1, F2[eFn].m_pdata,
                            g_n, g_np1_k, dt, mat, elasticity, solver_info);

      err_g = sqrt((g_np1_kp1 - g_np1_k)*(g_np1_kp1 - g_np1_k));
      g_np1_k = g_np1_kp1;      
    }
    err += Newton_Rapson4M(M.m_pdata, lambda, 
                           pFn.m_pdata, F2[pFnI].m_pdata, F2[Fa].m_pdata,F2[N].m_pdata,F2[A].m_pdata, 
                           g_np1_k, dt, mat, elasticity, solver_info,&is_it_cnvg_M,
                           &norm_R_0,d_gamma,k);                           

    if(err>0)
      break;
    
    // compute eFnp1 = Fr*M
    Matrix_AxB_3by3(F2+eFnp1,F2+Fa,&M); // <= Matrix_AxB(F2[eFnp1], 1.0,0.0,F2[Fa],0,M,0);
    err += Newton_Rapson_hardening(&g_np1_kp1, F2[eFnp1].m_pdata, 
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
    Matrix_inv(M,F2[MI]);
    Matrix_Tns2_AxBxC(pFnp1,1.0,0.0,F2[MI],pFn,F2[N]);  
    *g_out = g_np1_kp1;      
    *is_it_cnvg = 1;
  }
  else
    err += 1;
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);  
  free(F2); 
       
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

  enum {dF,F_s,Fn_s,pFn_s,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_init(double, F2[a],DIM_3,DIM_3,0.0);
    Matrix_eye(F2[a],DIM_3);
  }    
      
  double dt_s = dt/stepno;

  for(int b=0; b<DIM_3x3; b++)
  {
    F2[   dF].m_pdata[b] = (Fnp1_in[b] - Fn_in[b])/stepno;
    F2[ Fn_s].m_pdata[b] =  Fn_in[b];      
    F2[pFn_s].m_pdata[b] = pFn_in[b];
  }
  
  double g_n_s = g_n;

  int is_sub_cnvg = 0;        

  for(int a=0; a<stepno; a++)
  {
    if(DEBUG_PRINT_STAT)
      printf("-----------------------\n sub step: (%d/%d) \n-----------------------\n", a+1, stepno);
    
    is_sub_cnvg = 0;
    Matrix_AplusB(F2[F_s],1.0,F2[Fn_s],1.0,F2[dF]);

    err += staggered_Newton_Rapson_compute(pFnp1_out, M_out, g_out,lambda, 
                                         F2[pFn_s].m_pdata,
                                         F2[ Fn_s].m_pdata,
                                         F2[  F_s].m_pdata,
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
      F2[pFn_s].m_pdata[b] = pFnp1_out[b];
      F2[ Fn_s].m_pdata[b] = F2[F_s].m_pdata[b];      
    }
    g_n_s = *g_out;
  }
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);  

  *is_it_cnvg = is_sub_cnvg;

  free(F2);
    
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
  
  Matrix(double) AI, A, Fnp1, Fn;
  A.m_row    = A.m_col    = DIM_3;    A.m_pdata = I;
  Fnp1.m_row = Fnp1.m_col = DIM_3; Fnp1.m_pdata = Fnp1_out;
  Fn.m_row   = Fn.m_col   = DIM_3;   Fn.m_pdata = Fn_in;      

  Matrix_construct_redim(double,AI,DIM_3,DIM_3);
  Matrix_inv(A,AI);
  Matrix_AxB_3by3(&Fnp1,&AI,&Fn); // <= Matrix_AxB(Fnp1,1.0,0.0,AI,0,Fn,0);

  Matrix_cleanup(AI);  
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
