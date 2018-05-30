

#include"constitutive_model.h"
#include"poro_visco_plasticity.h"

constexpr const int Err = 1;

Tensor<2> I  = ttl::identity(i,j);
      
template<int R, int D = 10, class S = double>
using Matrix_10x10  = ttl::Tensor<R, D, S>;

double Macaulay(const double x,
                const double x0){
  return (x >= x0)? x : x0;  
}

/// class carrying 2nd order tensor (3x3 matrix) from C
class F_FromC
{
  public:
    double *Fnp1, *Fn, *pFnp1, *pFn, *hFnp1, *hFn;
    double *pc_np1;  
    
    F_FromC(){
      Fnp1 = Fn = pFnp1 = pFn = hFnp1 = hFn =NULL;
    }    
};
    

class StateVariables
{
  public:
    Tensor<2> M, MI, eFnp1, tau, eS, hat_eS, hat_tau, psi_d, psi_v;;
    Tensor<4> L;
    double pJ, pi, bar_tau, g_tau, g_pi, pi_m;
    
    StateVariables(){};
};

class PvpElasticity
{
  public:

    template<class T1, class T2, class T3> void compute_dWdC_dev(T1 &dWdC,
                                                                 T2 &d2WdC2,
                                                                 T3 &C,
                                                                 const double mu,
                                                                 bool compute_2nd_order = false){

      double CJ = ttl::det(C); // eC = eJ*eJ; pow(eCJ, -1.0/3.0) = pow(eJ, -2.0/3.0)
      double factor = pow(CJ, -1.0/3.0);
      
      double trC = C(i,i);
      
      Tensor<2> CI;
      int err = inv(C, CI);
      
      dWdC(i,j) = 0.5*mu*factor*(I(i,j) - 1.0/3.0*trC*CI(i,j));
    
      if(compute_2nd_order)
      {
        Tensor<4> CIxCI, dCIdC;
        CIxCI(i,j,k,l) = CI(i,j)*CI(k,l);
        dCIdC(i,j,k,l)  = -CI(i,k)*CI(j,l);
        
        d2WdC2(i,j,k,l) = 0.25*mu*factor*(-2.0/3.0*CI(i,j)*I(k,l) - 2.0/3.0*I(i,j)*CI(k,l) 
                                         + 2.0/9.0*trC*CIxCI(i,j,k,l) - 2.0/3.0*trC*dCIdC(i,j,k,l));
      }
      if(err>0)
        throw Err;
        
    }
    
    // U = 1/2*K*(J-1)*ln(J) + c*ln(J) + (alpha1 + alpha2*ln(J))*exp(-beta*ln(J))
    double compute_dUdJ(const double J,
                        const double K,
                        const double c,
                        const double alpha1,
                        const double alpha2,
                        const double beta){
      return 0.5*K*(log(J) + (J-1)/J) + c/J - ((alpha1 + alpha2*log(J))*beta - alpha2)*exp(-beta*log(J))/J;
    }
    
    double compute_d2UdJ2(const double J,
                          const double K,
                          const double c,
                          const double alpha1,
                          const double alpha2,
                          const double beta){
      return 0.5*K*(1.0/J + 1.0/J/J) - c/J/J + ((1.0+beta)*beta*(alpha1 + alpha2*log(J))
                                                 -(1.0+2.0*beta)*alpha2)*exp(-beta*log(J))/J/J;
    }
    
    template<class T1, class T2, class T3> void compute_elasticity_tensor(T1 &eS,
                                                                          T2 &L,
                                                                          T3 &eF,
                                                                          const double mu,
                                                                          const double K,
                                                                          const double c,
                                                                          const double alpha1,
                                                                          const double alpha2,
                                                                          const double beta,
                                                                          bool compute_2nd_order = false){
      int err = 0;
      
      Tensor<2> dWdC, eC, eCI;
      eC(i,j) = eF(k,i)*eF(k,j);
      
      err += inv(eC, eCI);
      double eJ = ttl::det(eF);
            
      compute_dWdC_dev(dWdC, L, eC, mu, compute_2nd_order);
        
      
      double dUdJ = compute_dUdJ(eJ, K, c, alpha1, alpha2, beta);
      
      eS(i,j) = 2.0*dWdC(i,j) + eJ*dUdJ*eCI(i,j);

      if(compute_2nd_order){
        double d2UdJ2 = compute_d2UdJ2(eJ, K, c, alpha1, alpha2, beta);
        L(i,j,k,l) = 4.0*L(i,j,k,l) + (eJ*dUdJ + eJ*eJ*d2UdJ2)*eCI(i,j)*eCI(k,l) - 2.0*eJ*dUdJ*eCI(i,k)*eCI(j,l);
      }
      
      if(err>0)
        throw Err;
    }
};

class PvpIntegrator
{  
  public:
    const MaterialPoroViscoPlasticity *mat;
    
    F_FromC Fs;    
    Tensor<2> Fr, FrI, eFn, Fa, N, A;
    
    StateVariables sv;
    PvpElasticity elasticity;
  
    double pc_np1;
    double pc_n;
    double dt;
    int max_itr;    
    double tol;
    double computer_zero;    
            
    PvpIntegrator(){
    }
    
    void set_pvp_material(const MaterialPoroViscoPlasticity *mat_in){mat=mat_in;};
    
    // set tensors, no memory is allocated
    void set_tenosrs(double *Fnp1,
                     double *Fn,
                     double *pFnp1,
                     double *pFn){      
      Fs.Fnp1  = Fnp1;
      Fs.Fn    = Fn;
      Fs.pFnp1 = pFnp1;
      Fs.pFn   = pFn;
      Fs.hFnp1 = Fs.hFn = I.data;
      compute_tensors();
    }
    
    // set tensors, no memory is allocated
    void set_tenosrs(double *Fnp1,
                     double *Fn,
                     double *pFnp1,
                     double *pFn,
                     double *hFnp1,
                     double *hFn){
      Fs.Fnp1  = Fnp1;
      Fs.Fn    = Fn;
      Fs.pFnp1 = pFnp1;
      Fs.pFn   = pFn;
      Fs.hFnp1 = hFnp1;
      Fs.hFn   = hFn;
      compute_tensors();
    }
    
    void set_numerical_parameters(const int max_itr_in = 6,
                                  const double tol_in = 1.0e-6,
                                  const double computer_zero_in = 1.0e-15){                                  
      computer_zero = computer_zero_in;
      max_itr = max_itr_in;
      tol     = tol_in;
    }
    void compute_tensors(void)
    {
      
      int err = 0;
      TensorA<2> pFn(Fs.pFn), hFn(Fs.hFn), Fn(Fs.Fn), Fnp1(Fs.Fnp1);
 
      Tensor<2> pFnI, hFnI, FnI;      
      err += inv(pFn, pFnI);
      err += inv(hFn, hFnI);
      err += inv(Fn,  FnI);

      // compute eFn
      eFn(i,j) = Fn(i,k)*hFnI(k,l)*pFnI(l,j);
      Fr(i,j)  = Fnp1(i,k)*FnI(k,j);
      Fa(i,j)  = Fr(i,k)*eFn(k,j);

      err += inv(Fr, FrI);

      // compute N : N = hFn*hFnp1I
      TensorA<2> hFnp1(Fs.hFnp1);
      Tensor<2>  hFnp1I;
      err += inv(hFnp1, hFnp1I);
      N(i,j) = hFn(i,k)*hFnp1I(k,j);
 
      // compute A : A = pFn*N_I*pFn_I
      Tensor<2>  NI;
            
      err += inv(N, NI);
      A(i,j) = pFn(i,k)*NI(k,l)*pFnI(l,j);
      
      Tensor<2> eFnI;
      err += inv(eFn, eFnI);
      sv.M(i,j) = eFnI(i,k)*FrI(k,l)*eFn(l,j);
            
      if(err>0)
        throw Err;        
    }
    
    
    // set scalars
    void set_scalars(double pc_np1_in,
                     double pc_n_in,
                     double dt_in){
      pc_np1 = pc_np1_in;
      pc_n   = pc_n_in;
      dt     = dt_in;
    }
    
    void update_StateVariables(const double pc){
      
      int err = 0;

      err += inv(sv.M, sv.MI);
      TensorA<2> pFnp1(Fs.pFnp1), pFn(Fs.pFn);
      
      pFnp1(i,j) = sv.MI(i,k)*pFn(k,l)*N(l,j);
      sv.eFnp1(i,j) = Fa(i,k)*sv.MI(k,j);
      
      sv.pJ = ttl::det(pFnp1);
      
      double c = compute_c(pc);
      double d = compute_d(pc);
      double K = (mat->K_p0 + c)*(d-1.0/d)/mat->K_kappa;
      double mu = mat->mu_0 + c*(d - 1.0/d)/mat->mu_1;
      double coeff_U_alpha = (mat->K_p0 + c)*pow(mat->c_inf/(mat->c_inf - c), mat->pl_n)*mat->K_kappa;
      double coeff_U_beta  = pow((mat->c_inf - c)/mat->c_inf, mat->pl_n)/mat->K_kappa;

      elasticity.compute_elasticity_tensor(sv.eS, sv.L, sv.eFnp1, mu, K, c, coeff_U_alpha, 0.0, coeff_U_beta, true);
      
      double tr_eS = sv.eS(i,i);
      sv.hat_eS(i,j) = sv.eS(i,j) - 1.0/3.0*tr_eS*I(i,j);
      sv.tau(i,j)    = sv.pJ*sv.eFnp1(i,k)*sv.eS(k,l)*sv.eFnp1(j,l);
      sv.pi          = compute_pi(sv.tau);
      
      double tr_tau = sv.tau(i,i);
      sv.hat_tau(i,j)  = sv.tau(i,j) - 1.0/3.0*tr_tau*I(i,j);
      sv.bar_tau = sqrt(sv.hat_tau(i,j)*sv.hat_tau(i,j));
      
      double a = compute_a(pc, c);
      sv.pi_m  = compute_pi_m(a, c);
      double b = compute_b(sv.pi, sv.pi_m);
      
      sv.g_tau = compute_g_tau(pc, a);
      sv.g_pi  = compute_g_pi(pc, a, b);
      
      compute_psi_d(sv.psi_d);
      compute_psi_v(sv.psi_v, pc);            
    }
    
    double compute_a(const double pc,
                     const double c){
      return (pc + c)/(1.0 + mat->yf_alpha);
    }
    
    double compute_c(const double pc){
      return mat->c_inf*(1.0-exp(-mat->c_Gamma*Macaulay(pc, mat->d_pcb)));
    }
    
    double compute_b(const double pi,
                     const double pi_m){
      return (pi >= pi_m)? mat->yf_alpha : 1;
    }
    
    double compute_d(const double pc){
      return 1.0 + mat->d_B*Macaulay(pc, mat->d_pcb);
    }
          
    double compute_g_tau(const double pc,
                         const double a){
      return sqrt(3.0/2.0)*a*mat->yf_M;
    }
    
    double compute_g_pi(const double pc,
                        const double a,
                        const double b){
      return a*b;
    }
    
    double compute_pi_m(const double a,
                        const double c){
      return a - c;
    }
    
    template<class T> double compute_pi(const T &tau){
      return -tau(i,i)/3.0;
    }
    
    double compute_H(const double pc)
    {
      return mat->hr_a1*exp(-mat->hr_Lambda1/pc) + mat->hr_a2*exp(-mat->hr_Lambda2/pc);
    }
    double compute_dHdp(const double pc){
      return 1.0/pc/pc*(mat->hr_a1*mat->hr_Lambda1*exp(-mat->hr_Lambda1/pc) + mat->hr_a2*mat->hr_Lambda2*exp(-mat->hr_Lambda2/pc));
    }
    
    double compute_beta_d(void){return 0.0;}
    double compute_beta_v(const double pc){
      return mat->cf_g0*(1.0 - pc/mat->cf_pcinf);
    }
    
    template<class T> void compute_psi_d(T &psi_d){

      double m_hat_eS = sqrt(sv.hat_eS(i,j)*sv.hat_eS(i,j));
      double beta_d = compute_beta_d();
      
      if(m_hat_eS<computer_zero)
        psi_d(i,j) = 1.0/3.0*beta_d*I(i,j);
      else
        psi_d(i,j) = sv.hat_eS(i,j)/m_hat_eS + 1.0/3.0*beta_d*I(i,j);
    }
    
    double compute_gamma_dot_d(const double d){
      return mat->flr_gamma_dot_0*(1.0-1.0/d)*pow(sv.bar_tau/sv.g_tau, 1.0/mat->flr_m);
    }
    double compute_gamma_dot_v(void){
      return (sv.pi>sv.pi_m)? mat->flr_gamma_dot_0*pow((sv.pi-sv.pi_m)/sv.g_pi, 1.0/mat->flr_m): 0.0;
    }
    
    template<class T> void compute_psi_v(T &psi_v,
                                         const double pc){
      double beta_v = compute_beta_v(pc);
      psi_v(i,j) = 1.0/3.0*beta_v*I(i,j);
    }
    /// compute_residual RM
    template<class T> void compute_RM(T &RM,
                                      const double pc){
      double d = compute_d(pc);
      double gamma_dot_d = compute_gamma_dot_d(d);                                        
      double gamma_dot_v = compute_gamma_dot_v();
      Tensor<2> pD;
      pD(i,j) = gamma_dot_d*sv.psi_d(i,j) + gamma_dot_v*sv.psi_v(i,j);
      
      RM(i,j) = pD(i,j) - 1.0/dt*(I(i,j) - sv.M(i,j));
    }
    
    /// compute residual Rpc
    /// returns Rpc
    void compute_Rp(double &Rpc,
                    const double pc){
      TensorA<2> pFn(Fs.pFn);
      Tensor<2> MIpFnN;
 
      MIpFnN(i,j) = sv.MI(i,k)*pFn(k,l)*N(l,j);
      double JM = ttl::det(MIpFnN);
      Rpc = log(JM) - compute_H(pc);
    }
    
    template<class T> void compute_d_bar_tau_d_tau(T &d_bar_tau){

      if(sv.bar_tau<computer_zero)
        d_bar_tau = {};
      else
        d_bar_tau(i,j) = sv.hat_tau(i,j)/sv.bar_tau;      
    }
    
    template<class T> void compute_dSdM(T &dSdM){
      Tensor<2> FaTFa;
      
      FaTFa(i,j) = Fa(k,i)*Fa(k,j);                                  
      dSdM(i,j,k,l) = 0.5*sv.L(i,j,m,r)*(I(m,l)*FaTFa(k,v)*sv.M(v,r) + sv.M(s,m)*FaTFa(s,k)*I(r,l));
    }
    template<class T1, class T2> void compute_d_tau_dM(T1 &d_tau_dM,
                                                       const T2 &dSdM){
      Tensor<4> deF_dM, deFT_dM;     

      deF_dM(i,j,k,l)  = Fa(i,k)*I(j,l);
      deFT_dM(i,j,k,l) = I(i,l)*Fa(j,k);
                  
      d_tau_dM(i,j,k,l) = -sv.tau(i,j)*sv.MI(l,k) + sv.pJ*deF_dM(i,m,k,l)*sv.eS(m,r)*sv.eFnp1(j,r)
                                                  + sv.pJ*sv.eFnp1(i,m)*dSdM(m,r,k,l)*sv.eFnp1(j,r)
                                                  + sv.eFnp1(i,m)*sv.eS(m,r)*deFT_dM(j,r,k,l);
    }
    
    template<class T> void compute_dgamma_dot_d_dM(T &dgamma_dot_d_dM,
                                                   const Tensor<4> &d_tau_dM,
                                                   const double d){
      double factor = mat->flr_gamma_dot_0*(1.0 - 1.0/d)/mat->flr_m*pow(sv.bar_tau/sv.g_tau, 1.0/mat->flr_m-1.0)/sv.g_tau;
      Tensor<2> d_bar_tau;
      compute_d_bar_tau_d_tau(d_bar_tau);
      
      dgamma_dot_d_dM(k,l) = factor*d_bar_tau(i,j)*d_tau_dM(i,j,k,l);
    }
    
    template<class T> void compute_dgamma_dot_v_dM(T &dgamma_dot_v_dM,
                                                   const Tensor<4> &d_tau_dM){
      double factor = 0;
      if(sv.pi>sv.pi_m)                                              
        factor = -1.0/3.0*mat->flr_gamma_dot_0/mat->flr_m/sv.g_pi*pow((sv.pi-sv.pi_m)/sv.g_pi, 1.0/mat->flr_m - 1.0);
      dgamma_dot_v_dM(k,l) = factor*I(i,j)*d_tau_dM(i,j,k,l);
    }
    
    template<class T1, class T2> void compute_dpsi_d_dM(T1 &dpsi_d_dM,
                                                        const T2 &dSdM){
      double m_hat_eS = sqrt(sv.hat_eS(i,j)*sv.hat_eS(i,j));
      
      if(m_hat_eS<computer_zero)
        dpsi_d_dM = {};
      else
        dpsi_d_dM(i,j,k,l) = 1.0/m_hat_eS*(dSdM(i,j,k,l) - I(i,j)*I(m, r)*dSdM(m,r,k,l))
                            -1.0/m_hat_eS/m_hat_eS/m_hat_eS*sv.hat_eS(i,j)*sv.hat_eS(m,r)*dSdM(m,r,k,l);
    }
          
    template<class T1, class T2> void compute_dRMdM(T1 &dRMdM,
                                                    T2 &dSdM,
                                                    const double pc){
      double d = compute_d(pc);
      Tensor<2> dgamma_dot_d_dM, dgamma_dot_v_dM;
      Tensor<4> d_tau_dM, dpsi_d_dM;
      
      compute_d_tau_dM(d_tau_dM, dSdM);
      
      double gamma_dot_d = compute_gamma_dot_d(d);
      
      compute_dgamma_dot_d_dM(dgamma_dot_d_dM, d_tau_dM,d);
      compute_dgamma_dot_v_dM(dgamma_dot_v_dM, d_tau_dM);

      // compute dpsi_d_dM
      compute_dpsi_d_dM(dpsi_d_dM, dSdM);
      
      dRMdM(i,j,k,l) = 1.0/dt*I(i,k)*I(j,l) + sv.psi_d(i,j)*dgamma_dot_d_dM(k,l) + gamma_dot_d*dpsi_d_dM(i,j,k,l) 
                                            + sv.psi_v(i,j)*dgamma_dot_v_dM(k,l);
    }
   
    double compute_dcdp(const double pc){
      return mat->c_inf*exp(-mat->c_Gamma*Macaulay(pc, mat->d_pcb));
    }
    
    double compute_dddp(const double pc){
      return (pc>mat->d_pcb)? mat->d_B: 0.0;
    }
    
    double compute_dadpc(const double pc){
      return 1.0/(1.0+mat->yf_alpha)*(1.0 + mat->cf_pcinf*exp(-mat->c_Gamma*Macaulay(pc, mat->d_pcb)));
    }
    
    double compute_dg_tau_dp(const double pc){
      return sqrt(3.0/2.0)*compute_dadpc(pc)*mat->yf_M;
    }
    
    double compute_dpJdp(const double pc){
      return -compute_dHdp(pc)*sv.pJ;
    }
    
    template<class T> void compute_dSdp(T &dSdp,
                                        const double pc){
      double c = compute_c(pc);
      double d = compute_d(pc);
            
      double dcdp = compute_dcdp(pc);
      double ddpc = compute_dddp(pc);
      double dmudp = (dcdp*(d-1.0/d) + c*(ddpc + 1.0/d/d*ddpc))*mat->mu_1;
      double dKdp  = (dcdp*(d-1.0/d) + (mat->K_p0 + c)*(ddpc + 1.0/d/d*ddpc))/mat->K_kappa;
      
      double pow_of_c = pow(mat->c_inf/(mat->c_inf - c), mat->pl_n);

      double coeff_U_alpha1 = dcdp*pow_of_c*mat->K_kappa*(1.0 + mat->pl_n*(mat->K_p0 + c)/(mat->c_inf - c));
      double coeff_U_alpha2 = dcdp*(mat->K_p0 + c)/(mat->c_inf - c)*mat->pl_n;
            
      double coeff_U_beta  = pow((mat->c_inf - c)/mat->c_inf, mat->pl_n)/mat->K_kappa;

      Tensor<4> L;
      elasticity.compute_elasticity_tensor(dSdp, L , sv.eFnp1, dmudp, dKdp, dcdp, coeff_U_alpha1, coeff_U_alpha2, coeff_U_beta, false);
    }
    
    template<class T1, class T2> void compute_d_tau_dp(T1 &d_tau_dp,
                                                       const T2 &dSdp,
                                                       const double pc){
      double dpJdp = compute_dpJdp(pc);            
      d_tau_dp(i,j) = sv.eFnp1(i,k)*(dpJdp*sv.eS(k,l) + sv.pJ*dSdp(k,l))*sv.eFnp1(j,l);
    }
    
    template<class T> double compute_d_bar_tau_dp(T &d_tau_dp){
      Tensor<2> d_bar_tau;
      
      compute_d_bar_tau_d_tau(d_bar_tau);
      return d_bar_tau(i,j)*d_tau_dp(i,j);
      
    }
    
    template<class T> double compute_d_bar_tau_g_tau_dp(const T &d_tau_dp,
                                                        const double pc){
      double d_bar_tau_dp = compute_d_bar_tau_dp(d_tau_dp);
      
      double d_g_tau_dp = compute_dg_tau_dp(pc);      
      return 1.0/sv.g_tau*d_bar_tau_dp - 1.0*sv.bar_tau/sv.g_tau/sv.g_tau*d_g_tau_dp;
    }
    
    template<class T> void compute_dRMdp(T &dRMdp,
                                         const double pc){
      // compute dSdp d_tau_dp
      double d = compute_d(pc);
      
      Tensor<2> dSdp, d_tau_dp;;
      compute_dSdp(dSdp, pc);
      compute_d_tau_dp(d_tau_dp, dSdp, pc);
      
      double bar_tau_g_tau = sv.bar_tau/sv.g_tau;
      double pow_bar_tau_g_tau = pow(sv.bar_tau/sv.g_tau, 1.0/mat->flr_m);      
                                                
      // compute dgamma_dot_d_pc
      double factor1 = mat->flr_gamma_dot_0/d/d*pow_bar_tau_g_tau;
      
      double dddp = compute_dddp(pc);      
      double d_bar_tau_g_tau_dp = compute_d_bar_tau_g_tau_dp(dSdp,pc);

      double dgamma_dot_d_pc = factor1*dddp;      
      if(sv.bar_tau>computer_zero)        
      {
        double factor2 = mat->flr_gamma_dot_0*(1.0-1.0/d)/mat->flr_m*pow_bar_tau_g_tau/bar_tau_g_tau;
        dgamma_dot_d_pc += factor2*d_bar_tau_g_tau_dp;
      }
      
      // compute gamma_dot_d
      double gamma_dot_d = compute_gamma_dot_d(d);
      
      // compute d_psi_d_dp
      Tensor<2> d_psi_d_dp = {};
      double m_hat_eS = sqrt(sv.hat_eS(i,j)*sv.hat_eS(i,j));

      if(m_hat_eS>computer_zero){        
        double tr_dSdp = dSdp(i,i);      
        double hat_eS_dSdp = 1.0/m_hat_eS/m_hat_eS/m_hat_eS*sv.hat_eS(i,j)*dSdp(i,j);
        d_psi_d_dp(i,j) = 1.0/m_hat_eS*(dSdp(i,j) - tr_dSdp*I(i,j)) - hat_eS_dSdp*sv.hat_eS(i,j);
      }
      
      // compute dgamma_dot_v_pc
      double factor3 = 0.0;
      if(sv.pi>sv.pi_m)
        factor3 = mat->flr_gamma_dot_0/mat->flr_m*pow((sv.pi-sv.pi_m)/sv.g_pi, 1.0/mat->flr_m-1.0);

      double dpidp = -1.0/3.0*I(i,j)*d_tau_dp(i,j);
      double dg_pidp = compute_b(sv.pi, sv.pi_m)*compute_dadpc(pc);
      double dgamma_dot_v_pc = factor3*(dpidp/sv.g_pi - (sv.pi-sv.pi_m)/sv.g_pi/sv.g_pi*dg_pidp);

            
      // compute gamma_dot_v
      double gamma_dot_v = compute_gamma_dot_v();
      
      // compute d_psi_v_dp
      Tensor<2> d_psi_v_dp;
      d_psi_v_dp = 1.0/3.0*mat->cf_g0/mat->cf_pcinf*I(i,j);
            
      
      dRMdp(i,j) = dgamma_dot_d_pc*sv.psi_d(i,j) + gamma_dot_d*d_psi_d_dp(i,j) + dgamma_dot_v_pc*sv.psi_v(i,j) + gamma_dot_v*d_psi_v_dp(i,j);
      
    }
    
    void compute_left_matrix(Matrix_10x10<2> &K,
                             const double pc){
      Tensor<4> dRMdM, dSdM;
      Tensor<2> dRMdp, dRpdM;
                  
      compute_dSdM(dSdM);
      compute_dRMdM(dRMdM, dSdM, pc);                    
      compute_dRMdp(dRMdp, pc);

      dRpdM(i,j) = -sv.MI(j,i);      
      double dRpdp = -compute_dHdp(pc);
      
      for(int ia=0; ia<DIM_3x3; ia++){
        K(ia,DIM_3x3) = dRMdp.get(ia);
        K(DIM_3x3,ia) = dRpdM.get(ia);
        for(int ib=0; ib<DIM_3x3; ib++){
          K(ia,ib) = dRMdM.get(ia*DIM_3x3+ib);
        }
      }
      K(DIM_3x3,DIM_3x3) = dRpdp;      
    }
    
    void compute_right_vector(Matrix_10x10<1> &R,
                              const double pc){
      
      Tensor<2> RM;
      compute_RM(RM, pc);
      double Rpc = 0.0;
      compute_Rp(Rpc, pc);
      for(int ia=0; ia<DIM_3x3; ia++)
        R.data[ia] = RM.data[ia];
      
      R.data[DIM_3x3] = Rpc;
    }
    
    ~PvpIntegrator(){};    
};

int poro_visco_plasticity_integration_algorithm(const MaterialPoroViscoPlasticity *mat,
                                                double *Fnp1,
                                                double *Fn,
                                                double *pFnp1,
                                                double *pFn,
                                                double *pc_np1,
                                                const double pc_n,
                                                const double dt)
{
  int err = 0;
  
  PvpIntegrator pvp;
  pvp.set_pvp_material(mat);
  pvp.set_tenosrs(Fnp1, Fn, pFnp1, pFn);
  pvp.set_scalars(*pc_np1, pc_n, dt);
  pvp.set_numerical_parameters(1000);
  
  Matrix_10x10<1> R, du;
  Matrix_10x10<2> K, KI;
  
  Tensor<2> dM;
  
  
  double pc = pc_n;
  double dpc = 0.0;
  
  double norm_R_0 = pvp.computer_zero;
        
  for(int iA = 0; iA<pvp.max_itr; iA++){
    
    if(iA==0){
      pvp.update_StateVariables(pc);  
      pvp.compute_right_vector(R, pc);
      norm_R_0 = sqrt(R(i)*R(i));
    }
  
    pvp.compute_left_matrix(K, pc);
  
    try
    {
      KI(i,j) = ttl::inverse(K)(i,j);  // attempt to take the inverse
    }
    catch(const int inverseException)  // no inverse exists
    {
      if(1)
        printf( "Matrix is singular. The solution (Poro-visco-plasticity) could not be computed.\n");

      ++err;
      break;
    } 

//  print_tensor(K, "K");
//  print_tensor(KI, "KI");
  
    du(i) = -KI(i,j)*R(j);
  
    for(int ia=0; ia<DIM_3; ++ia)
      for(int ib=0; ib<DIM_3; ++ib)
        dM.data[ia*DIM_3 + ib] = du.get(ia*DIM_3 + ib);

    dpc = du.get(DIM_3x3);
    
    pvp.sv.M(i,j) = pvp.sv.M(i,j) + dM(i,j);
    if(dpc<0)
    {
      printf("dpc = %e\n", dpc);  
      dpc = 0.0;
    }  
    pc = pc + dpc;
    
    pvp.update_StateVariables(pc);  
    pvp.compute_right_vector(R, pc);    
    double norm_R = sqrt(R(i)*R(i));
  
    printf("(%d/%d)residual: |R| = %e, |R0| = %e, |R/R0| = %e\n", iA, pvp.max_itr, norm_R, norm_R_0, norm_R/norm_R_0);
    if(norm_R/norm_R_0<pvp.tol)
      break;
  } 
  *pc_np1 = pc;
  
  return err; 
}