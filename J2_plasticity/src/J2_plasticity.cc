#include "constitutive_model.h"
#include "constitutive_model_handle.h"   //for exascale variables
#include "J2_plasticity.h"
#include "hyperelasticity.h"

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81
#define J2P_INT_ALG_TOL 1.0e-10

#define idx_2(row,col) ((unsigned int) ((row)*3+(col)))

template <class T, class U, class V>
static int compute_push_forward(const T& F,
                                const U& A,
                                V& a)
{
  int err = 0;
  a = F(i,k)*A(k,l)*F(j,l);
  return err;
}

template <class T, class U, class V>
static int compute_pull_back(const T& F,
                             const U& a,
                             V& A)
{
  int err = 0;
  Tensor<2> FI;
  err += inv(F, FI);
  compute_push_forward(FI, a, A);
  return err;
}

template <class T, class U, class V>
static int compute_pull_Tensor4(const T& FI,
                                const U& aep,
                                V& Aep)
{
  int err = 0;
  Aep(i,j,k,l) = FI(i,m)*FI(j,r)*FI(k,v)*FI(l,w)*aep(m,r,v,w);
  return err;
}

// dev(a) = a - 1/3 tr(a) i
template <class T, class U>
static int compute_dev(const T& a,
                       U& dev_a)
{
  int err = 0;
  double tra = a(i,i)/3.0;
  Tensor<2> eye = ttl::identity(i,j);

  dev_a(i,j) = a(i,j) - tra*eye(i,j);;
  return err;
}

template <class T, class U>
static int compute_s0(double G,
                      const T& bbar,
                      U& s0)
{
  int err = 0;
  double tra = bbar(i,i)/3.0;
  Tensor<2> eye = ttl::identity(i,j);

  s0(i,j) = G*bbar(i,j) - G*tra*eye(i,j);
  return err;
}

/* bbar = J^-(2/3) F F' */
template <class T, class U>
static int compute_bbar_of_J(T& bbar,
                             const U& F,
                             double J)
{
  int err = 0;
  double J23 = pow(J,-2.0/3.0);
  bbar(i,j) = J23*F(i,k)*F(j,k);
  return err;
}

template <class T, class U>
static int compute_bbar(T& bbar,
                        const U& F)
{
  int err = 0;

  double J = ttl::det(F);
  err += compute_bbar_of_J(bbar,F,J);
  return err;
}

template <class T, class U, class V>
static int compute_Fubar(const T& F,
                         const U& Fn,
                         V& Fubar)
{
  int err = 0;
  Tensor<2> FnI;
  err += inv(Fn, FnI);
  Fubar(i,j) = F(i,k)*FnI(k,j);
  double J = ttl::det(Fubar);
  double Ju13 = pow(J, -1.0/3.0);

  for (int ia = 0; ia < DIM_3x3; ia++) {
    Fubar.data[ia] *= Ju13;
  }

  return err;
}

template <class T, class U, class V, class W>
static int compute_sp_tr(const T& F,
                         const U& Fn,
                         const V& spn,
                         W& sp_tr)
{
  int err = 0;
  Tensor<2> Fubar;

  err += compute_Fubar(F,Fn,Fubar);
  err += compute_push_forward(Fubar,spn,sp_tr);

  double tr = sp_tr(i,i)/3.0;
  sp_tr[0][0] -= tr;
  sp_tr[1][1] -= tr;
  sp_tr[2][2] -= tr;

  return err;
}

template <class T, class U, class V>
static double compute_normal(const MATERIAL_J2_PLASTICITY *J2P,
                             const MATERIAL_ELASTICITY *mat_e,
                             const T& s_tr,
                             const U& sp_tr,
                             V& n)
{
  double G = mat_e->G;
  double coef = (J2P->hp)/(3.0*G)*(1.0 - J2P->beta);

  /* compute ksi_tr and ||ksi_tr|| simultaneously */
  n(i,j) = s_tr(i,j) - coef*sp_tr(i,j);
  double nrm = n(i,j)*n(i,j);
  nrm = sqrt(nrm);

  /* compute normal: n = ksi_tr / ||ksi_tr|| */
  if(nrm > 0)
    for(int ia = 0; ia < DIM_3x3; ia++) n.data[ia] /= nrm;

  return nrm;
};

template <class T, class U, class V>
static double phi_yield_functioncompute_normal(const MATERIAL_J2_PLASTICITY *J2P,
                                               const MATERIAL_ELASTICITY *mat_e,
                                               const T& s_tr,
                                               const U& sp_tr,
                                               V& n,
                                               double ep_n)
{
  double ksi_nrm = compute_normal(J2P,mat_e,s_tr,sp_tr,n);
  return ksi_nrm - sqrt(2.0/3.0)*(J2P->k0 + (J2P->beta)*(J2P->hp)*ep_n);
}

int J2_plasticity_integration_alg(double *sp_out,
                                  double *ep_out,
                                  double *gamma_out,
                                  double *F_in,
                                  double *Fn_in,
                                  double *sp_n_in,
                                  double ep_n,
                                  MATERIAL_J2_PLASTICITY *J2P,
                                  MATERIAL_ELASTICITY *mat_e)
{
  int err = 0;

  perIter_ODE_EXA_metric += 2;  //2 ODEs per J2 plasticity
  double G = mat_e->G;

  TensorA<2> F(F_in), Fn(Fn_in), sp(sp_out), sp_n(sp_n_in);
  Tensor<2>  bbar, s_tr, sp_tr, n, s0;

  double J = ttl::det(F);
  double J23 = pow(J,-2.0/3.0);

  // compute bbar at n + 1
  err += compute_bbar_of_J(bbar, F,J);

  // compute mu_bar
  double mu_bar = G*J23*bbar(i,i)/3.0;

  // compute sp_tr
  err += compute_sp_tr(F,Fn,sp_n,sp_tr);

  // compute s_tr
  err += compute_s0(G, bbar, s0);
  s_tr(i,j) = s0(i,j) - sp_tr(i,j);

  // compute ksi_tr and the normal of plastic loading

  double phi = phi_yield_functioncompute_normal(J2P,mat_e,s_tr,sp_tr,n, ep_n);

  if(phi <= J2P_INT_ALG_TOL)
  {
    *gamma_out = 0.0;
    *ep_out = ep_n;
    sp(i,j) = sp_tr(i,j);
  }
  else
  {
    *gamma_out = phi/(2.0* mu_bar*(1.0 + J2P->hp/(3.0*G)*(1.0 - J2P->beta)
                            + J2P->beta*J2P->hp/(3.0*mu_bar)));
    double tmp = 2.0*mu_bar*(*gamma_out);
    *ep_out = ep_n + sqrt(2.0/3.0) *(*gamma_out);

    sp(i,j) = sp_tr(i,j) + tmp*n(i,j);
  }

  return err;
}

template <class T, class U, class V, class W>
static int compute_S0_Sbar(T& S0,
                           U& Sbar,
                           const V& F,
                           const W& sp,
                           HyperElasticity *elast)
{
  int err = 0;
  // compute the current configuration deviatoric stresses

  double G = (elast->mat)->G;
  double kappa = (elast->mat)->kappa;

  Tensor<2> s0,sbar,bbar, C, CI;

  err += compute_bbar(bbar, F);
  err += compute_s0(G, bbar,s0);
  sbar(i,j) = s0(i,j) - sp(i,j);

  //perform pull-back */
  err += compute_pull_back(F,s0,S0);
  err += compute_pull_back(F,sbar,Sbar);

  //compute volumetric stress
  C(i,j) = F(k,i)*F(k,j);
  err += inv(C,CI);
  double J = ttl::det(F);
  double dudj = 0.0;
  elast->compute_dudj(&dudj,J);

  S0(i,j) = S0(i,j) + kappa*J*dudj*CI(i,j);
  Sbar(i,j) = Sbar(i,j) + kappa*J*dudj*CI(i,j);

  return err;
}

template <class T, class U, class V, class W>
static int compute_S0_Sbar_dev(T& S0,
                               U& Sbar,
                               const V& F,
                               const W& sp,
                               HyperElasticity *elast)
{
  int err = 0;
  // compute the current configuration deviatoric stresses

  double G = (elast->mat)->G;
  Tensor<2> s0,sbar,bbar, C, CI;

  err += compute_bbar(bbar, F);
  err += compute_s0(G, bbar,s0);
  sbar(i,j) = s0(i,j) - sp(i,j);

  //perform pull-back */
  err += compute_pull_back(F,s0,S0);
  err += compute_pull_back(F,sbar,Sbar);
  return err;
}

int compute_S0_Sbar_public(double *S0_out,
                           double *Sbar_out,
                           double *F_in,
                           double *sp_in,
                           HyperElasticity *elast)
{
  int err = 0;

  TensorA<2> S0(S0_out), Sbar(Sbar_out), F(F_in), sp(sp_in);
  err += compute_S0_Sbar(S0,Sbar,F,sp,elast);

  return err;
}

int compute_S0_Sbar_dev_public(double *S0_out,
                               double *Sbar_out,
                               double *F_in,
                               double *sp_in,
                               HyperElasticity *elast)
{
  int err = 0;

  TensorA<2> S0(S0_out), Sbar(Sbar_out), F(F_in), sp(sp_in);
  err += compute_S0_Sbar_dev(S0,Sbar,F,sp,elast);

  return err;
}

int compute_S0_Sbar_split_public(double *dS0_out,   double *vS0_out,
                                 double *dSbar_out, double *vSbar_out,
                                 double *F_in, double *sp_in,
                                 HyperElasticity *elast)
{
  int err = 0;

  TensorA<2> dS0(dS0_out), vS0(vS0_out), dSbar(dSbar_out), vSbar(vSbar_out), F(F_in), sp(sp_in);

  // compute the current configuration deviatoric stresses
  double G = (elast->mat)->G;
  double kappa = (elast->mat)->kappa;

  Tensor<2> s0,sbar,bbar, C, CI;

  err += compute_bbar(bbar, F);
  err += compute_s0(G, bbar,s0);
  sbar(i,j) = s0(i,j) - sp(i,j);

  //perform pull-back */
  err += compute_pull_back(F,s0,dS0);
  err += compute_pull_back(F,sbar,dSbar);

  //compute volumetric stress
  C(i,j) = F(k,i)*F(k,j);
  err += inv(C,CI);
  double J = ttl::det(F);
  double dudj = 0.0;
  elast->compute_dudj(&dudj,J);
  double tmp = kappa*J*dudj;
  vS0(i,j)   = tmp*CI(i,j);
  vSbar(i,j) = tmp*CI(i,j);

  return err;
}

// compute the deviatoric initial/unloading tangent in the reference
// configuration
template <class T, class U, class V, class W>
static int compute_unloading_Aep_dev(T& Aep_dev,
                                     const U& F,
                                     const V& Fn,
                                     const W& sp_n,
                                     double G)
{
  int err = 0;
  Tensor<2> C, CI, Spn, eye = ttl::identity(i,j);


  // compute C and related terms
  //compute volumetric stress
  C(i,j) = F(k,i)*F(k,j);
  err += inv(C,CI);
  double J23 = ttl::det(C);
  J23 = pow(J23, -1.0/3.0);

  double Cpp = C(i,i);

  // compute pull-back of spn
  err += compute_pull_back(F,sp_n,Spn);

  double CSp = C(i,j)*Spn(i,j);
  Aep_dev(i,j,k,l) = 2.0/3.0*J23*(((G*Cpp-CSp)*(CI(i,k)*CI(j,l) + CI(i,j)*CI(k,l)/3.0))
                 - (CI(i,j)*(G*eye(k,l) - Spn(k,l)) + (G*eye(i,j) - Spn(i,j))*CI(k,l)));

  return err;
}


// compute the deviatoric plastic loading tangent in the reference
// configuration
template <class T, class U, class V, class W>
static int compute_loading_Aep_dev(T& Aep_dev,
                                   const U& F,
                                   const V& Fn,
                                   const W& sp_n,
                                   double gamma,
                                   const MATERIAL_J2_PLASTICITY *J2P,
                                   const MATERIAL_ELASTICITY *mat_e)
{
  int err = 0;
  double G = mat_e->G;
  Tensor<2> C,CI,Spn,bbar,devbbar,FI,sp_tr,Fubar,s_tr,normal,normal2,zeros;
  Tensor<2> eye = ttl::identity(i,j);

  C(i,j) = F(k,i)*F(k,j);
  err += inv(C,CI);
  double J23 = ttl::det(C);
  J23 = pow(J23, -1.0/3.0);

  err += inv(F,FI);
  err += compute_bbar(bbar, F);
  err += compute_dev(bbar, devbbar);

  // compute Itr = G tr(bbar) - tr(Fubar spn Fubar')
  double Itr = G*bbar(i,i);
  {
    err += compute_Fubar(F,Fn,Fubar);
    err += compute_push_forward(Fubar,sp_n,sp_tr);
    double trace_sp = sp_tr(i,i);
    sp_tr[0][0] -= trace_sp/3.0;
    sp_tr[1][1] -= trace_sp/3.0;
    sp_tr[2][2] -= trace_sp/3.0;
    Itr -= trace_sp;
  }

  // compute s_tr, normal and ||s_tr||
  err += compute_s0(G, bbar,s_tr);
  s_tr(i,j) = s_tr(i,j) - sp_tr(i,j);
  double norm_s_tr = s_tr(i,j)*s_tr(i,j);
  norm_s_tr = sqrt(norm_s_tr);

  compute_normal(J2P,mat_e,s_tr,sp_tr,normal);
  normal2(i,j) = normal(i,k)*normal(k,j);

  // compute factors
  double mu_bar = G*J23*bbar(i,i)/3.0;
  double f0 = 1.0 - 2.0*mu_bar*gamma/norm_s_tr;
  double del0 = 1.0 + J2P->hp/(3.0*mu_bar);
  double f1 = (1.0/del0 - 1.0 + f0);
  double del1 = f1 * 2.0/3.0*Itr;
  double del2 = (1.0/del0 - 1.0)*4.0/3.0*mu_bar*gamma - 2.0/3.0*norm_s_tr*f1;
  double del3 = 2.0*norm_s_tr*f1;
  double del4 = (1.0/del0 - 1.0)*4.0/3.0*G*gamma*J23;

  Tensor<4> aep = (f0*(2.0/3.0*Itr*(eye(i,k)*eye(l,j) - eye(i,j)*eye(k,l)/3.0)
                    - 2.0/3.0*(s_tr(i,j)*eye(k,l) + eye(i,j)*s_tr(k,l)))
                    - del1*normal(i,j)*normal(k,l)
                    - del2*(normal(i,j)*eye(k,l)     + normal(k,l)*eye(i,j))*0.5
                    - del3*(normal(i,j)*normal2(k,l) + normal(k,l)*normal2(i,j))*0.5
                    - del4*(normal(i,j)*devbbar(k,l) + normal(k,l)*devbbar(i,j))*0.5);

  err += compute_pull_Tensor4(FI,aep,Aep_dev);
  return err;
}

template <class T, class U, class V, class W>
static int compute_Lbar(T& Lbar,
                        const U& F,
                        const V& Fn,
                        const W& sp_n,
                        double gamma,
                        const MATERIAL_J2_PLASTICITY *J2P,
                        HyperElasticity *elast)
{
  int err = 0;
  double kappa = (elast->mat)->kappa;
  double G = (elast->mat)->G;
  double J = ttl::det(F);

  Tensor<2> C, CI;
  C(i,j) = F(k,i)*F(k,j);
  err += inv(C,CI);

  double dudj = 0.0;
  double d2udj2 = 0.0;

  /* compute deviatoric tangent */
  if (gamma > 0)
    err += compute_loading_Aep_dev(Lbar,F,Fn,sp_n,gamma,J2P,elast->mat);
  else
    err += compute_unloading_Aep_dev(Lbar,F,Fn,sp_n,G);

  elast->compute_dudj(&dudj,J);
  elast->compute_d2udj2(&d2udj2,J);

  double coeff_1 = kappa*J*(dudj+J*d2udj2);
  double coeff_2 = 2.0*kappa*J*dudj;

  Lbar(i,j,k,l) = ((coeff_1*CI(i,j)*CI(k,l))-(coeff_2*CI(i,k)*CI(l,j)));

  return err;
}


int compute_Lbar_public(double *Lbar_out,
                        double *F_in,
                        double *Fn_in,
                        double *sp_n_in,
                        double gamma,
                        MATERIAL_J2_PLASTICITY *J2P,
                        HyperElasticity *elast)
{
  int err = 0;

  TensorA<4> Lbar(Lbar_out);
  TensorA<2> F(F_in), Fn(Fn_in), sp_n(sp_n_in);

  err += compute_Lbar(Lbar, F, Fn, sp_n, gamma, J2P, elast);
  return err;
}

int compute_Lbar_dev_public(double *Lbar_out,
                            double *F_in,
                            double *Fn_in,
                            double *sp_n_in,
                            double gamma,
                            MATERIAL_J2_PLASTICITY *J2P,
                            HyperElasticity *elast)
{
  int err = 0;

  TensorA<4> Lbar(Lbar_out);
  TensorA<2> F(F_in), Fn(Fn_in), sp_n(sp_n_in);

  if (gamma > 0)
    err += compute_loading_Aep_dev(Lbar,F,Fn,sp_n,gamma,J2P,elast->mat);
  else
    err += compute_unloading_Aep_dev(Lbar,F,Fn,sp_n,(elast->mat)->G);
  return err;
}

int compute_Lbar_split_public(double *dLbar_out,
                              double *vLbar_out,
                              double *F_in,
                              double *Fn_in,
                              double *sp_n_in,
                              double gamma,
                              MATERIAL_J2_PLASTICITY *J2P,
                              HyperElasticity *elast)
{
  int err = 0;

  TensorA<4> dLbar(dLbar_out), vLbar(vLbar_out);
  TensorA<2> F(F_in), Fn(Fn_in), sp_n(sp_n_in);

  double kappa = (elast->mat)->kappa;
  double G = (elast->mat)->G;
  double J = ttl::det(F);

  Tensor<2> C, CI;
  C(i,j) = F(k,i)*F(k,j);
  err += inv(C,CI);

  double dudj = 0.0;
  double d2udj2 = 0.0;

  // compute deviatoric stiffness
  if (gamma > 0)
    err += compute_loading_Aep_dev(dLbar,F,Fn,sp_n,gamma,J2P,elast->mat);
  else
    err += compute_unloading_Aep_dev(dLbar,F,Fn,sp_n,G);

  elast->compute_dudj(&dudj,J);
  elast->compute_d2udj2(&d2udj2,J);

  double coeff_1 = kappa*J*(dudj+J*d2udj2);
  double coeff_2 = 2.0*kappa*J*dudj;

  vLbar(i,j,k,l) = ((coeff_1*CI(i,j)*CI(k,l))-(coeff_2*CI(i,k)*CI(l,j)));

  return err;
}

int J2_plasticity_update_elasticity(MATERIAL_J2_PLASTICITY *J2P,
                                    HyperElasticity *elast,
                                    double *F_in,
                                    double *Fn_in,
                                    double *sp_in,
                                    double *sp_n_in,
                                    double gamma,
                                    const int compute_stiffness)
{
  int err = 0;
  TensorA<2> F(F_in),Fn(Fn_in),sp(sp_in),sp_n(sp_n_in), S(elast->S);
  Tensor<2> S0;

  err += compute_S0_Sbar(S0,S,F,sp,elast);

  if(compute_stiffness)
  {
    TensorA<4> L(elast->L);
    err += compute_Lbar(L,F,Fn,sp_n,gamma,J2P,elast); //compute stiffness
  }

  return err;
}
