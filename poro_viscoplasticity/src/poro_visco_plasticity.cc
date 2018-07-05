/**
 * Authors:
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include"constitutive_model.h"
#include"poro_visco_plasticity.h"
#include"GcmSolverInfo.h"

constexpr const int Err = 1;

Tensor<2> I  = ttl::identity(i,j);
Tensor<2> O  = ttl::identity(i,j)*0.0;
      
template<int R, int D = 10, class S = double>
using Matrix_10x10  = ttl::Tensor<R, D, S>;

double Macaulay(const double x,
                const double x0){
  return (x >= x0)? x-x0 : 0.0;  
}

template<class T> void print(double *A, const int dim, const char *name = NULL)
{
  if(name==NULL)
    printf("A={");
  else
    printf("%s = {", name);
    
  for(int ia=0; ia<dim; ia++){
    if(ia==dim-1)
      printf("%.17e,", A[ia]);
    else
      printf("%.17e", A[ia]);
  }
  printf("};\n");
}

template<class T> void print(T &A, const int dim, const char *name = NULL)
{
  if(name==NULL)
    printf("A={");
  else
    printf("%s = {", name);
    
  for(int ia=0; ia<dim; ia++){
    if(ia==dim-1)
      printf("%.17e", A.data[ia]);
    else
      printf("%.17e, ", A.data[ia]);
  }
  printf("};\n");
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
    
class PvpMaterial
{
  public:
    const MaterialPoroViscoPlasticity *param;
    
    PvpMaterial(){param=NULL;};
    PvpMaterial(const MaterialPoroViscoPlasticity *param_in){set_pvp_material_parameters(param_in);};
    void set_pvp_material_parameters(const MaterialPoroViscoPlasticity *param_in){param=param_in;};
    
    /// compute shear modulus as a function of pc
    ///
    /// \param[in] pc
    /// \return    shear modulus
    double compute_shear_modulus(const double pc)
    {
      double c = compute_c(pc);
      double d = compute_d(pc);
      return compute_shear_modulus(c, d);      
    }

    /// compute shear modulus as a function of c d
    /// whenever c and d are available, this funciton allows to skip re computing c and d.
    ///
    /// \param[in] c pvp model function value of pc
    /// \param[in] d pvp model function value of pc
    /// \return    shear modulus 
    double compute_shear_modulus(const double c,
                                 const double d){
      return param->mu_0 + c*(d - 1.0/d)*param->mu_1;
    }
    
    /// compute bulk mudulus as a function of pc
    ///
    /// \param[in] pc
    /// \return    bulk modulus    
    double compute_bulk_modulus(const double pc){
      double c = compute_c(pc);
      double d = compute_d(pc);
      return compute_bulk_modulus(c, d);
    }

    /// compute bulk modulus as a function of c d
    /// whenever c and d are available, this funciton allows to skip re computing c and d.
    ///
    /// \param[in] c pvp model function value of pc
    /// \param[in] d pvp model function value of pc
    /// \return    bulk modulus        
    double compute_bulk_modulus(const double c,
                                const double d){
      return (param->K_p0 + c)*(d-1.0/d)/param->K_kappa;
    }

    /// compute material constant for volumetric part of strain energy density function
    /// U = 1/2*K*(J-1)*ln(J) + c*ln(J) + (alpha1 + alpha2*ln(J))*exp(-beta*ln(J))
    /// note: alpha2 = 0 for computing stress and elasticity tensor
    ///
    /// \param[in] c pvp model function value of pc
    /// \return    material constant (alpha_1)
    double compute_coeff_U_alpha(const double c){
      return (param->K_p0 + c)*pow(param->c_inf/(param->c_inf - c), param->pl_n)*param->K_kappa;
    }

    /// compute material constant for volumetric part of strain energy density function
    /// U = 1/2*K*(J-1)*ln(J) + c*ln(J) + (alpha1 + alpha2*ln(J))*exp(-beta*ln(J))
    /// note: alpha2 = 0 for computing stress and elasticity tensor
    ///
    /// \param[in] c pvp model function value of pc
    /// \return    material constant (beta)        
    double compute_coeff_U_beta(const double c){
      return pow((param->c_inf - c)/param->c_inf, param->pl_n)/param->K_kappa;
    } 
    
    /// compute a
    ///
    /// \param[in] pc
    /// \param[in] c pvp model function value of pc
    /// \return    a
    double compute_a(const double pc,
                     const double c){
      return (pc + c)/(1.0 + param->yf_alpha);
    }
    
    /// compute c
    ///
    /// \param[in] pc
    /// \return    c
    double compute_c(const double pc){
      return param->c_inf*(1.0-exp(-param->c_Gamma*Macaulay(pc, param->d_pcb)));
    }

    /// compute b
    ///
    /// \param[in] pi   Kirchhoff pressure(-tr(tau)/3.0)
    /// \param[in] pi_m parameter for yeild surface 
    /// \param[in] c    pvp model function value of pc
    /// \return    b   
    double compute_b(const double pi,
                     const double pi_m){
      return (pi >= pi_m)? param->yf_alpha : 1.0;
    }

    /// compute d
    ///
    /// \param[in] pc
    /// \return    d    
    double compute_d(const double pc){
      return 1.0 + param->d_B*Macaulay(pc, param->d_pcb);
    }

    /// compute g_tau
    ///
    /// \param[in] pc
    /// \param[in] a pvp model function value of pc
    /// \return    g_tau          
    double compute_g_tau(const double pc,
                         const double a){
      return sqrt(1.5)*a*param->yf_M; // sqrt(3.0/2.0) = sqrt(1.5)
    }

    /// compute g_pi
    ///
    /// \param[in] a pvp model function value of pc
    /// \param[in] b pvp model function value of pc
    /// \return    g_pi
    double compute_g_pi(const double a,
                        const double b){
      return a*b;
    }

    /// compute pi_m
    ///
    /// \param[in] a pvp model function value of pc
    /// \param[in] c pvp model function value of pc
    /// \return    pi_m    
    double compute_pi_m(const double a,
                        const double c){
      return a - c;
    }

    /// compute H = -a1*exp(-L1/pc) - a2*exp(-L2/pc);
    ///
    /// \param[in] pc
    /// \return    H     
    double compute_H(const double pc)
    {
      return -param->hr_a1*exp(-param->hr_Lambda1/pc) - param->hr_a2*exp(-param->hr_Lambda2/pc);
    }

    /// compute gamma_dot_d = gamma_dot_0*(1.0-1.0/d)*(bar_tau/g_tau)^(1.0/m);
    ///
    /// \param[in] d           pvp model function value of pc
    /// \param[in] bar_tau     tau_h/||tau_hat||
    /// \param[in] g_tau       parameter for yeild surface 
    /// \return    gamma_dot_d
    double compute_gamma_dot_d(const double d,
                               const double bar_tau,
                               const double g_tau){
      return param->flr_gamma_dot_0*(1.0-1.0/d)*pow(bar_tau/g_tau, 1.0/param->flr_m);
    }

    /// compute gamma_dot_v = gamma_dot_0*((pi-pi_m)/g_pi)^(1.0/m)
    ///                                 or
    ///                       0.0 if pi<=pi_m
    ///
    /// \param[in] pi          Kirchhoff pressure(-tr(tau)/3.0)
    /// \param[in] pi_m        parameter for yeild surface 
    /// \param[in] g_pi        parameter for yeild surface 
    /// \return    gamma_dot_v    
    double compute_gamma_dot_v(const double pi,
                               const double pi_m,
                               const double g_pi){
      return (pi>pi_m)? param->flr_gamma_dot_0*pow((pi-pi_m)/g_pi, 1.0/param->flr_m): 0.0;
    }

    /// compute dilatancy function value. zero for this model 
    /// \return 0.0 always 
    double compute_beta_D(void){return 0.0;}

    /// compute compaction function value 
    /// 
    /// \param[in] pc
    /// \return double     
    double compute_beta_C(const double pc){
      return param->cf_g0*(1.0 - pc/param->cf_pcinf);
    }        
};

/// state variable carrier
class StateVariables
{
  public:
    Tensor<2> M, MI, eFnp1, tau, eS, hat_eS, hat_tau, psi_d, psi_v;;
    Tensor<4> L;
    double pJ, pi, bar_tau, g_tau, g_pi, pi_m;
    
    StateVariables(){};
};

/// object for the pvp model integration
class PvpElasticity
{
  public:

    template<class T1, class T2, class T3> void compute_dWdC_dev(T1 &dWdC,
                                                                 T2 &d2WdC2,
                                                                 T3 &C,
                                                                 const double mu,
                                                                 bool compute_4th_order = false){

      double CJ = ttl::det(C); // eCJ = eJ*eJ; pow(eCJ, -1.0/3.0) = pow(eJ, -2.0/3.0)
      double factor = pow(CJ, -1.0/3.0);
      
      double trC = C(i,i);
      
      Tensor<2> CI;
      int err = inv(C, CI);
      
      dWdC(i,j) = 0.5*mu*factor*(I(i,j) - 1.0/3.0*trC*CI(i,j));
    
      if(compute_4th_order)
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
      return 0.5*K*(log(J) + (J-1)/J) + c/J + (alpha2 - beta*(alpha1 + alpha2*log(J)))*exp(-beta*log(J))/J;
    }
    
    double compute_d2UdJ2(const double J,
                          const double K,
                          const double c,
                          const double alpha1,
                          const double alpha2,
                          const double beta){
      double JJ = 1.0/J/J;
      return 0.5*K*(1.0/J + 1.0*JJ) - c*JJ + ((1.0+beta)*beta*(alpha1 + alpha2*log(J))
                                                 -(1.0+2.0*beta)*alpha2)*exp(-beta*log(J))*JJ;
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
                                                                          bool compute_4th_order = false){
      int err = 0;
      
      Tensor<2> dWdC, eC, eCI;
      eC(i,j) = eF(k,i)*eF(k,j);
      
      err += inv(eC, eCI);
      double eJ = ttl::det(eF);
            
      compute_dWdC_dev(dWdC, L, eC, mu, compute_4th_order);
        
      
      double dUdJ = compute_dUdJ(eJ, K, c, alpha1, alpha2, beta);
      
      eS(i,j) = 2.0*dWdC(i,j) + eJ*dUdJ*eCI(i,j);

      if(compute_4th_order){
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
    PvpMaterial    mat;
    StateVariables sv;
    PvpElasticity  elasticity;
        
    F_FromC Fs;    
    Tensor<2> Fr, FrI, eFn, Fa, N, A;
  
    double pc_np1;
    double pc_n;
    const GcmSolverInfo *solver_info;   
            
    PvpIntegrator(){
    }
    
    /// set material object
    ///
    /// \param[in] param_in a poniter of MaterialPoroViscoPlasticity object
    void set_pvp_material_parameters(const MaterialPoroViscoPlasticity *param_in){mat.set_pvp_material_parameters(param_in);};
    
    /// Set tensors. pointers are used to be able to support c. No memory is allocated.
    /// 
    /// \param[in] *Fnp1  total deformation gradient at t(n+1)
    /// \param[in] *Fn    total deformation gradient at t(n)
    /// \param[in] *pFnp1 plastic deformation gradient at t(n+1)
    /// \param[in] *pFn   plastic deformation gradient at t(n)
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
    
    /// Set tensors. pointers are used to be able to support c. No memory is allocated.
    /// 
    /// \param[in] *Fnp1  total deformation gradient at t(n+1)
    /// \param[in] *Fn    total deformation gradient at t(n)
    /// \param[in] *pFnp1 plastic deformation gradient at t(n+1)
    /// \param[in] *pFn   plastic deformation gradient at t(n)
    /// \param[in] *hFnp1 thermal expansion part of deformation gradient at t(n+1)
    /// \param[in] *hFn   thermal expansion part of deformation gradient at t(n)    
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
    
    /// senumerical parameters such as tolerance and maximum number of iterations
    ///
    /// \param[in] max_itr_in       maximum number of NR iterations
    /// \param[in] tol_in           NR tolerance
    /// \param[in] computer_zero_in computer zero
    void set_solver_info(const GcmSolverInfo *solver_info_in){
      solver_info = solver_info_in;
    }
    
    /// compute initial tensors such as eFn, Fr, Fa, N, Fa, A, and M
    /// if inverse has error, error exception will be thrown
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
      Tensor<2> NI;
            
      err += inv(N, NI);
      A(i,j) = pFn(i,k)*NI(k,l)*pFnI(l,j);
      
      Tensor<2> eFnI;
      err += inv(eFn, eFnI);
      sv.M(i,j) = eFnI(i,k)*FrI(k,l)*eFn(l,j);
            
      if(err>0)
        throw Err;        
    }
    
    
    /// set scalars pc and dt
    ///
    /// \param[in] pc_np1_in pc at t(n+1)
    /// \param[in] pc_n_in   pc at t(n)
    void set_scalars(double pc_np1_in,
                     double pc_n_in){
      pc_np1 = pc_np1_in;
      pc_n   = pc_n_in;
    }
    
    /// After new values(pFnp1, and pc) are computed, internal state variables can be updated 
    /// using this function. First, pc dependent values will be computed and stress, elasticity tensor, and
    /// tau, tau_bar will be updated accordingly. Computing sequence is ordered by their dependency.
    ///
    /// \param[in] pc
    void update_StateVariables(const double pc){
      
      int err = 0;

      err += inv(sv.M, sv.MI);
      TensorA<2> pFn(Fs.pFn), pFnp1(Fs.pFnp1);
      
      pFnp1(i,j) = sv.MI(i,k)*pFn(k,l)*N(l,j);

      sv.eFnp1(i,j) = Fa(i,k)*sv.M(k,j);
      sv.pJ = ttl::det(pFnp1);
      
      double c  = mat.compute_c(pc);
      double d  = mat.compute_d(pc);
      double K  = mat.compute_bulk_modulus(c, d);
      double mu = mat.compute_shear_modulus(c, d);

      double coeff_U_alpha = mat.compute_coeff_U_alpha(c);
      double coeff_U_beta  = mat.compute_coeff_U_beta(c);

      elasticity.compute_elasticity_tensor(sv.eS, sv.L, sv.eFnp1, mu, K, c, coeff_U_alpha, 0.0, coeff_U_beta, true);
      
      double tr_eS = sv.eS(i,i);
      sv.hat_eS(i,j) = sv.eS(i,j) - 1.0/3.0*tr_eS*I(i,j);
      sv.tau(i,j)    = sv.pJ*sv.eFnp1(i,k)*sv.eS(k,l)*sv.eFnp1(j,l);
      sv.pi          = compute_pi(sv.tau);
      
      double tr_tau = sv.tau(i,i);
      sv.hat_tau(i,j)  = sv.tau(i,j) - 1.0/3.0*tr_tau*I(i,j);
      sv.bar_tau = sqrt(sv.hat_tau(i,j)*sv.hat_tau(i,j));
      
      double a = mat.compute_a(pc, c);
      sv.pi_m  = mat.compute_pi_m(a, c);
      double b = mat.compute_b(sv.pi, sv.pi_m);
      
      sv.g_tau = mat.compute_g_tau(pc, a);
      sv.g_pi  = mat.compute_g_pi(a, b);
      
      compute_psi_d(sv.psi_d);
      compute_psi_v(sv.psi_v, pc);            
    }
    
    template<class T> double compute_pi(const T &tau){
      return -tau(i,i)/3.0;
    }
    
    double compute_dHdp(const double pc){
      return -1.0/(pc*pc)*(mat.param->hr_a1*mat.param->hr_Lambda1*exp(-mat.param->hr_Lambda1/pc) + mat.param->hr_a2*mat.param->hr_Lambda2*exp(-mat.param->hr_Lambda2/pc));
    }    
    
    template<class T> void compute_psi_d(T &psi_d){

      double m_hat_eS = sqrt(sv.hat_eS(i,j)*sv.hat_eS(i,j));
      double beta_D = mat.compute_beta_D();
      
      if(m_hat_eS<solver_info->computer_zero)
        psi_d(i,j) = 1.0/3.0*beta_D*I(i,j);
      else
        psi_d(i,j) = sv.hat_eS(i,j)/m_hat_eS + 1.0/3.0*beta_D*I(i,j);
    }
    
    template<class T> void compute_psi_v(T &psi_v,
                                         const double pc){
      double beta_C = mat.compute_beta_C(pc);
      psi_v(i,j) = -1.0/3.0*beta_C*I(i,j);
    }
    /// compute_residual RM
    template<class T> void compute_RM(T &RM,
                                      const double pc){
      Tensor<2> pD;
      compute_pD(pD, pc);      
      RM(i,j) = solver_info->dt*pD(i,j) - I(i,j) + A(i,k)*sv.M(k,j);
    }
    
    /// compute residual Rpc
    /// returns Rpc
    void compute_Rp(double &Rpc,
                    const double pc){
                      
      Rpc = sv.pJ - exp(mat.compute_H(pc));
                            
      //TensorA<2> pFn(Fs.pFn);
      //Tensor<2> MIpFnN;
 
      //MIpFnN(i,j) = sv.MI(i,k)*pFn(k,l)*N(l,j);
      //double JM = ttl::det(MIpFnN);
      //Rpc = log(JM) - mat.compute_H(pc);
    }
    
    template<class T> void compute_d_bar_tau_d_tau(T &d_bar_tau){

      if(sv.bar_tau<solver_info->computer_zero)
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
      //deFT_dM(i,j,k,l) = I(i,l)*Fa(j,k);
                  
      d_tau_dM(i,j,k,l) = -sv.tau(i,j)*sv.MI(l,k) + sv.pJ*deF_dM(i,m,k,l)*sv.eS(m,r)*sv.eFnp1(j,r)
                                                  + sv.pJ*sv.eFnp1(i,m)*dSdM(m,r,k,l)*sv.eFnp1(j,r)
                                                  + sv.pJ*sv.eFnp1(i,m)*sv.eS(m,r)*deF_dM(j,r,k,l);
    }
    
    template<class T> void compute_dgamma_dot_d_dM(T &dgamma_dot_d_dM,
                                                   const Tensor<4> &d_tau_dM,
                                                   const double d){
      double factor = mat.param->flr_gamma_dot_0/(mat.param->flr_m*sv.g_tau)*(1.0 - 1.0/d)*pow(sv.bar_tau/sv.g_tau, 1.0/mat.param->flr_m-1.0);
      Tensor<2> d_bar_tau;
      compute_d_bar_tau_d_tau(d_bar_tau);
      
      dgamma_dot_d_dM(k,l) = factor*d_bar_tau(i,j)*d_tau_dM(i,j,k,l);
    }
    
    template<class T> void compute_dgamma_dot_v_dM(T &dgamma_dot_v_dM,
                                                   const Tensor<4> &d_tau_dM){
      if(sv.pi>sv.pi_m){        
        double factor = -mat.param->flr_gamma_dot_0/(3.0*mat.param->flr_m*sv.g_pi)*pow((sv.pi-sv.pi_m)/sv.g_pi, 1.0/mat.param->flr_m - 1.0);
        dgamma_dot_v_dM(k,l) = factor*I(i,j)*d_tau_dM(i,j,k,l);
      }
      else
        dgamma_dot_v_dM(k,l) = O(k,l);

    }
    
    template<class T1, class T2> void compute_dpsi_d_dM(T1 &dpsi_d_dM,
                                                        const T2 &dSdM){
      double m_hat_eS = sqrt(sv.hat_eS(i,j)*sv.hat_eS(i,j));
      
      if(m_hat_eS<solver_info->computer_zero)
        dpsi_d_dM = {};
      else
        dpsi_d_dM(i,j,k,l) = 1.0/m_hat_eS*(dSdM(i,j,k,l) - I(i,j)*I(m, r)*dSdM(m,r,k,l)/3.0)
                            -1.0/m_hat_eS/m_hat_eS/m_hat_eS*sv.hat_eS(i,j)*sv.hat_eS(m,r)*dSdM(m,r,k,l);
    }
          
    template<class T1, class T2> void compute_dRMdM(T1 &dRMdM,
                                                    T2 &dSdM,
                                                    const double pc){
      double d = mat.compute_d(pc);
      Tensor<2> dgamma_dot_d_dM, dgamma_dot_v_dM;
      Tensor<4> d_tau_dM, dpsi_d_dM;
      
      compute_d_tau_dM(d_tau_dM, dSdM);
      
      double gamma_dot_d = mat.compute_gamma_dot_d(d, sv.bar_tau, sv.g_tau);
      
      compute_dgamma_dot_d_dM(dgamma_dot_d_dM, d_tau_dM,d);
      compute_dgamma_dot_v_dM(dgamma_dot_v_dM, d_tau_dM);

      // compute dpsi_d_dM
      compute_dpsi_d_dM(dpsi_d_dM, dSdM);
      
      dRMdM(i,j,k,l) = A(i,k)*I(j,l) + solver_info->dt*(sv.psi_d(i,j)*dgamma_dot_d_dM(k,l) + gamma_dot_d*dpsi_d_dM(i,j,k,l) 
                                            + sv.psi_v(i,j)*dgamma_dot_v_dM(k,l));
    }
   
    double compute_dcdp(const double pc){
      return mat.param->c_inf*mat.param->c_Gamma*exp(-mat.param->c_Gamma*Macaulay(pc, mat.param->d_pcb));
    }
    
    double compute_dddp(const double pc){
      return (pc>mat.param->d_pcb)? mat.param->d_B: 0.0;
    }
    
    double compute_dadpc(const double pc){
      return 1.0/(1.0+mat.param->yf_alpha)*(1.0 + compute_dcdp(pc));
    }
    
    double compute_dg_tau_dp(const double pc){
      return sqrt(1.5)*compute_dadpc(pc)*mat.param->yf_M; // sqrt(3.0/2.0) = sqrt(1.5)
    }
    
    template<class T> void compute_dSdp(T &dSdp,
                                        const double pc){
      double c = mat.compute_c(pc);
      double d = mat.compute_d(pc);
            
      double dcdp = compute_dcdp(pc);
      double ddpc = compute_dddp(pc);
      double dmudp = (dcdp*(d-1.0/d) + c*(1.0 + 1.0/(d*d))*ddpc)*mat.param->mu_1;
      double dKdp  = (dcdp*(d-1.0/d) + (mat.param->K_p0 + c)*(1.0 + 1.0/(d*d))*ddpc)/mat.param->K_kappa;
      
      double pow_of_c = pow(mat.param->c_inf/(mat.param->c_inf - c), mat.param->pl_n);

      double coeff_U_alpha1 = dcdp*pow_of_c*mat.param->K_kappa*(1.0 + mat.param->pl_n*(mat.param->K_p0 + c)/(mat.param->c_inf - c));
      double coeff_U_alpha2 = dcdp*(mat.param->K_p0 + c)/(mat.param->c_inf - c)*mat.param->pl_n;
            
      double coeff_U_beta  = mat.compute_coeff_U_beta(c);

      Tensor<4> L;
      elasticity.compute_elasticity_tensor(dSdp, L , sv.eFnp1, dmudp, dKdp, dcdp, coeff_U_alpha1, coeff_U_alpha2, coeff_U_beta, false);
    }
    
    template<class T1, class T2> void compute_d_tau_dp(T1 &d_tau_dp,
                                                       const T2 &dSdp,
                                                       const double pc){
      double dHdp = compute_dHdp(pc);                                 
//      double dpJdp = sv.pJ*dHdp(pc);          
//      d_tau_dp(i,j) = sv.eFnp1(i,k)*(dpJdp*sv.eS(k,l) + sv.pJ*dSdp(k,l))*sv.eFnp1(j,l);
      d_tau_dp(i,j) = sv.eFnp1(i,k)*sv.pJ*dSdp(k,l)*sv.eFnp1(j,l) - dHdp*sv.tau(i,j);
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
      return 1.0/sv.g_tau*d_bar_tau_dp - sv.bar_tau/sv.g_tau/sv.g_tau*d_g_tau_dp;
    }
    
    template<class T> void compute_dRMdp(T &dRMdp,
                                         const double pc){
      // compute dSdp d_tau_dp
      double d = mat.compute_d(pc);
      
      Tensor<2> dSdp, d_tau_dp;
      compute_dSdp(dSdp, pc);
      compute_d_tau_dp(d_tau_dp, dSdp, pc);
      
      double bar_tau_g_tau = sv.bar_tau/sv.g_tau;
      double pow_bar_tau_g_tau = pow(sv.bar_tau/sv.g_tau, 1.0/mat.param->flr_m);      
                                                
      // compute dgamma_dot_d_pc
      double factor1 = mat.param->flr_gamma_dot_0/(d*d)*pow_bar_tau_g_tau;
      
      double dddp = compute_dddp(pc);      
      double d_bar_tau_g_tau_dp = compute_d_bar_tau_g_tau_dp(d_tau_dp,pc);

      double dgamma_dot_d_pc = factor1*dddp;      
      if(sv.bar_tau>solver_info->computer_zero)        
      {
        double factor2 = mat.param->flr_gamma_dot_0*(1.0-1.0/d)/mat.param->flr_m*pow_bar_tau_g_tau/bar_tau_g_tau;
        dgamma_dot_d_pc += factor2*d_bar_tau_g_tau_dp;
      }
      
      // compute gamma_dot_d
      double gamma_dot_d = mat.compute_gamma_dot_d(d, sv.bar_tau, sv.g_tau);
      
      // compute d_psi_d_dp
      Tensor<2> d_psi_d_dp = {};
      double m_hat_eS = sqrt(sv.hat_eS(i,j)*sv.hat_eS(i,j));

      if(m_hat_eS>solver_info->computer_zero){        
        double tr_dSdp = dSdp(i,i);      
        double hat_eS_dSdp = 1.0/m_hat_eS/m_hat_eS/m_hat_eS*sv.hat_eS(i,j)*dSdp(i,j);
        d_psi_d_dp(i,j) = 1.0/m_hat_eS*(dSdp(i,j) - tr_dSdp*I(i,j)/3.0) - hat_eS_dSdp*sv.hat_eS(i,j);
      }
      
      // compute dgamma_dot_v_pc
      double factor3 = 0.0;
      if(sv.pi>sv.pi_m)
        factor3 = mat.param->flr_gamma_dot_0/(mat.param->flr_m*sv.g_pi)*pow((sv.pi-sv.pi_m)/sv.g_pi, 1.0/mat.param->flr_m-1.0);

      double dadpc = compute_dadpc(pc);
      double dcdp = compute_dcdp(pc);
      double dpidp = -d_tau_dp(i,i)/3.0;
      double dpi_mdp = dadpc - dcdp;
      double dg_pidp = mat.compute_b(sv.pi, sv.pi_m)*dadpc;
      double dgamma_dot_v_pc = factor3*(dpidp*0.0 - dpi_mdp - (sv.pi-sv.pi_m)/sv.g_pi*dg_pidp);

      // compute gamma_dot_v
      double gamma_dot_v = mat.compute_gamma_dot_v(sv.pi, sv.pi_m, sv.g_pi);
      
      // compute d_psi_v_dp
      Tensor<2> d_psi_v_dp = mat.param->cf_g0/(3.0*mat.param->cf_pcinf)*I(i,j);

      dRMdp(i,j) = solver_info->dt*(dgamma_dot_d_pc*sv.psi_d(i,j) + gamma_dot_d*d_psi_d_dp(i,j) + dgamma_dot_v_pc*sv.psi_v(i,j) + gamma_dot_v*d_psi_v_dp(i,j));      
    }
    
    template<class T> void compute_pD(T &pD,
                                      const double pc){
      double d = mat.compute_d(pc);
      double gamma_dot_d = mat.compute_gamma_dot_d(d, sv.bar_tau, sv.g_tau);
      double gamma_dot_v = mat.compute_gamma_dot_v(sv.pi, sv.pi_m, sv.g_pi);
      pD(i,j) = gamma_dot_d*sv.psi_d(i,j) + gamma_dot_v*sv.psi_v(i,j);                                          
    }
    
    void compute_left_matrix(Matrix_10x10<2> &K,
                             const double pc){
      Tensor<4> dRMdM, dSdM;
      Tensor<2> dRMdp, dRpdM;
                  
      compute_dSdM(dSdM);
      compute_dRMdM(dRMdM, dSdM, pc);                    
      compute_dRMdp(dRMdp, pc);

      dRpdM(i,j) = -sv.pJ*sv.MI(j,i);    

      double Hpc = mat.compute_H(pc);                
      double dRpdp = -sv.pJ*compute_dHdp(pc);
      
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
    
    double compute_pc(double pcn,
                      const double pJ){
      
      double logJp = log(pJ);      
      double DfDpcn = compute_dHdp(pcn);
      
      int it=0;
      double pcnp1 = pc_n;
      
      if(fabs(DfDpcn) < solver_info->tol_M){
        double a = pcn;
        double b = a;
        double c = a;
        
        while(1){
          b = b*2.0;
          double fpcnp1 = logJp + mat.compute_H(b);
          if(fpcnp1 <= 0)
             break;
        }
        
        const int maxIt=500;
        
        while (it < maxIt){
          ++it;
          c = 0.5 * (a+b);
          double fpccheck = -logJp + mat.compute_H(c);
          if(fabs(fpccheck) < solver_info->tol_hardening)
            break;
          else if( fpccheck > 0 )
            a = c;
          else
            b =c;
        }
        pcnp1 = c;
      } else {
        const int maxIt=100;
        pcnp1 = pcn;
        
        while(it < maxIt){
          ++it;
          double f = -logJp + mat.compute_H(pcnp1);
          double dHdp = compute_dHdp(pcnp1);
          
          double dpc = -f/dHdp;
          pcnp1 = pcnp1 + dpc;          
          if(fabs(logJp - mat.compute_H(pcnp1)) < solver_info->tol_hardening)
            break;
        }
      }
      return pcnp1;                        
    }

    template<class T> void compute_dMdu(double* dMdUs,
                                        const double pc,
                                        double* Grad_us,
                                        const int nne,
                                        const int ndofn){
      Tensor<4> U;
      Tensor<2> B, dSdp;
      Tensor<4> II, chi, eC, IoxI, eSdoxeSd;

      II(i,j,k,l) = I(i,k)*I(j,l);
      IoxI(i,j,k,l) = I(i,j)*I(k,l);
      eSdoxeSd(i,j,k,l) = sv.hat_eS(i,j)*sv.hat_eS(k,l);
      
      eC(i,j) = sv.eFnp1(k,i)*sv.eFnp1(k,j);
      
      double d = mat.compute_d(pc);
      double dddp = compute_dddp(pc);
      double dHdp = compute_dHdp(pc);
      double one_over_dH = 1.0/dHdp;
      
      // start computing A(left side);
      // compute chi
      Tensor<2> eFFa_sym;
      eFFa_sym(i,j) = 0.5*(sv.eFnp1(k,i)*Fa(k,j) + Fa(k,i)*sv.eFnp1(k,j));
      
      compute_dSdp(dSdp, pc);
      chi(i,j,k,l) = sv.L(i,j,k,m)*eFFa_sym(m,l) - one_over_dH*dSdp(i,j)*sv.MI(l,k);
            
      //compute D(gamma_dot_d) part      
      double gamma_dot_d = mat.compute_gamma_dot_d(d, sv.bar_tau, sv.g_tau);
      double gamma_dot_v = mat.compute_gamma_dot_v(sv.pi, sv.pi_m, sv.g_pi);      

      double bar_tau_g_tau = sv.bar_tau/sv.g_tau;
      double pow_bar_tau_g_tau = pow(sv.bar_tau/sv.g_tau, 1.0/mat.param->flr_m);
      double d_g_tau_dp = compute_dg_tau_dp(pc); 
      
      double factor1 = pow_bar_tau_g_tau*one_over_dH*mat.param->flr_gamma_dot_0*(-dddp/(d*d) 
                       + (1.0 - 1.0/d)/(mat.param->flr_m*sv.g_tau)*d_g_tau_dp);
      double factor2 = sv.pJ*mat.param->flr_gamma_dot_0/(mat.param->flr_m*sv.g_tau)*(1.0 - 1.0/d)*pow_bar_tau_g_tau/sv.bar_tau*sv.g_tau;

      Tensor<2> U1, U2, d_bar_tau_d_tau, eFeSeFT, d_bar_tau_d_tau_sym;
      Tensor<2> eFTdtaueF;
      eFTdtaueF(i,j) = sv.eFnp1(k,i)*d_bar_tau_d_tau(k,l)*sv.eFnp1(l,j);
            
      compute_d_bar_tau_d_tau(d_bar_tau_d_tau);
      d_bar_tau_d_tau_sym(i,j) = 0.5*(d_bar_tau_d_tau(i,j) + d_bar_tau_d_tau(j,i));
      
      eFeSeFT(i,j) = sv.eFnp1(i,k)*sv.eS(k,l)*sv.eFnp1(j,l);
      
      double sub_factor1 = d_bar_tau_d_tau(i,j)*eFeSeFT(i,j);

      U1(i,j) = factor1*sv.MI(j,i);
      U2(i,j) = factor2*(-sub_factor1*sv.MI(j,i) + 2.0*Fa(k,i)*d_bar_tau_d_tau_sym(k,l)*sv.eFnp1(l,m)*sv.eS(m,j)
                         + eFTdtaueF(k,l)*chi(k,l,i,j));
                         
      // compute D(psi_d) part
      Tensor<4> U3;      
      double m_eSd = sqrt(sv.hat_eS(i,j)*sv.hat_eS(i,j));
      double one_over_m_eSd = 0.0;
      if(m_eSd>solver_info->computer_zero)
      U3(i,j,k,l) = one_over_m_eSd*chi(i,j,k,l) + (IoxI(i,j,r,s) + eSdoxeSd(i,j,r,s))*chi(r,s,k,l);
      
      // compute D(gamma_dot_v) part
      double factor3 = mat.param->flr_gamma_dot_0/(mat.param->flr_m*sv.g_pi)*pow((sv.pi-sv.pi_m)/sv.g_pi, 1.0/mat.param->flr_m - 1.0);
      Tensor<2> U4, U5, U6;
      U4(i,j) = sv.pJ/3.0*I(k,l)*eFeSeFT(k,l)*sv.MI(j,i);
      U5(i,j) = -sv.pJ/3.0*(2.0*Fa(k,i)*sv.eFnp1(k,l)*sv.eS(l,j) + eC(k,l)*chi(k,l,i,j));

      double dadpc = compute_dadpc(pc);
      double dcdp = compute_dcdp(pc);      
      double dpi_mdp = dadpc - dcdp;
      double dg_pidp = mat.compute_b(sv.pi, sv.pi_m)*dadpc;
      
      U6(i,j) = one_over_dH*(dpi_mdp + (sv.pi-sv.pi_m)/sv.g_pi*dg_pidp)*sv.MI(j,i);
      
      // compute D(psi_v) part
      double factor4 = mat.param->cf_g0/3.0/mat.param->cf_pcinf*one_over_dH;
            
      U(i,j,k,l) = A(i,k)*I(j,l) + solver_info->dt*sv.psi_d(i,j)*(U1(k,l) + U2(k,l))
                                 + solver_info->dt*gamma_dot_d*U3(i,j,k,l)
                                 + solver_info->dt*factor3*sv.psi_v(i,j)*(U4(k,l) + U5(k,l) + U6(k,l))
                                 - solver_info->dt*gamma_dot_v*factor4*I(i,j)*sv.MI(l,k);
                                      
      // start computing B(right side);
      
      for(int ia=0; ia<nne; ia++){
        for(int ib=0; ib<ndofn; ib++){
          int id_ab = ia*ndofn*DIM_3x3 + DIM_3x3*ib;
          
          TensorA<2> dMdu(dMdUs + id_ab), Grad_u(Grad_us + id_ab);
          
          Tensor<2> GradeFnMeSeF, GradeFnMeSeF_sym, eFGradeFnM, eFGradeFnM_sym;
          GradeFnMeSeF(i,j) = Grad_u(i,k)*eFn(k,l)*sv.M(l,m)*sv.eFnp1(j,m);
          eFGradeFnM(i,j) = sv.eFnp1(k,i)*Grad_u(k,l)*eFn(l,m)*sv.M(m,j);
          
          GradeFnMeSeF_sym(i,j) = 0.5*(GradeFnMeSeF(i,j) + GradeFnMeSeF(j,i));
          eFGradeFnM_sym(i,j) = 0.5*(eFGradeFnM(i,j) + eFGradeFnM(j,i));
          
          Tensor<4> U7;
          U7(i,j,k,l) = one_over_m_eSd*sv.L(i,j,k,l) + (IoxI(i,j,r,s)+ eSdoxeSd(i,j,r,s))*sv.L(r,s,k,l);
          double factor5 = d_bar_tau_d_tau(i,j)*GradeFnMeSeF_sym(i,j) + eFTdtaueF(i,j)*sv.L(i,j,k,l)*eFGradeFnM_sym(k,l);
          double factor6 = I(i,j)*GradeFnMeSeF_sym(i,j) + eC(i,j)*sv.L(i,j,k,l)*eFGradeFnM_sym(k,l);
          B(i,j) = -solver_info->dt*factor2*factor5*sv.psi_d(i,j) - solver_info->dt*gamma_dot_v*U7(i,j,k,l)*eFGradeFnM_sym(k,l)
                   -solver_info->dt*sv.pJ/3.0*factor3*factor6*sv.psi_v(i,j);
          
          try{
            dMdu = ttl::solve( U, B);
          }
          catch(int solve_err){
            dMdu(i,j) = I(i,j);
            printf("error on computing dMdu for pvp model. dMdu = delta_ij is set. Solution may not be converging.\n");
          }     
          
        }
      }

    }

    ~PvpIntegrator(){};    
};

int poro_visco_plasticity_integration_algorithm_staggered(const MaterialPoroViscoPlasticity *mat,
                                                          const GcmSolverInfo *solver_info,
                                                          double *Fnp1_in,
                                                          double *Fn_in,
                                                          double *pFnp1_in,
                                                          double *pFn_in,
                                                          double *pc_np1,
                                                          const double pc_n){
  int err = 0;
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(mat);       
  pvp.set_tenosrs(Fnp1_in, Fn_in, pFnp1_in, pFn_in); 

  pvp.set_scalars(*pc_np1, pc_n);
  pvp.set_solver_info(solver_info);
                                                            
  TensorA<2> pFnp1(pFnp1_in);
  
  pvp.update_StateVariables(pc_n);    
  double pcr = pc_n;
      
  int it=0;
  
  // staggered solution for M      
  while (it < solver_info->max_itr_stag){
    it++;
    int itMr=0;
    while(itMr < solver_info->max_itr_stag){
      
     Tensor<2> rm;
     pvp.compute_RM(rm, pcr);
      
     double norm_rm = sqrt(rm(i,j)*rm(i,j));
     double norm_Mr = sqrt(pvp.sv.M(i,j)*pvp.sv.M(i,j));
      
      if ( norm_rm / norm_Mr < solver_info->tol_M  )
        break;
      
      // if convergence is not achieved, iterate
      
      itMr++;
      
      Tensor<4> dRMdM, dSdM;
                  
      pvp.compute_dSdM(dSdM);
      pvp.compute_dRMdM(dRMdM, dSdM, pcr);  
      
      Tensor<2> dm = {};
      
      try
      {      
        dm = ttl::solve(dRMdM, rm);      
      }
      catch (int i)
      {
        ++err;
        printf("Error detected\n");
        break;      
      }
            
      pvp.sv.M(i,j) = pvp.sv.M(i,j) - dm(i,j);
      pvp.update_StateVariables(pcr);        
    }
        
    double pJ = ttl::det(pFnp1);
    *pc_np1 = pvp.compute_pc(pc_n, pJ);
        
    pcr = *pc_np1;
    
    Tensor<2> rm;
    pvp.update_StateVariables(pcr);
    pvp.compute_RM(rm, pcr);
    
    double norm_rm = sqrt(rm(i,j)*rm(i,j));
    double norm_Mr = sqrt(pvp.sv.M(i,j)*pvp.sv.M(i,j));    
    
    if (norm_rm/norm_Mr < solver_info->tol_M  )
      break;    
  }
  return err;
}

int poro_visco_plasticity_integration_algorithm_implicit_00(const MaterialPoroViscoPlasticity *mat,
                                                            const GcmSolverInfo *solver_info,
                                                            double *Fnp1,
                                                            double *Fn,
                                                            double *pFnp1,
                                                            double *pFn,
                                                            double *pc_np1,
                                                            const double pc_n)
{
  int err = 0;
  
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(mat);
  pvp.set_tenosrs(Fnp1, Fn, pFnp1, pFn);
  pvp.set_scalars(*pc_np1, pc_n);
  pvp.set_solver_info(solver_info);
  
  Matrix_10x10<1> R, du;
  Matrix_10x10<2> K = {}, KI = {};
  
  Tensor<2> dM;
  
  double pc = pc_n;
  double dpc = 0.0;
  
  double norm_R_0 = pvp.solver_info->computer_zero;
  double norm_R   = pvp.solver_info->computer_zero;
  double eng_norm_0 = pvp.solver_info->computer_zero;
        
  for(int iA = 0; iA<pvp.solver_info->max_itr_M; iA++){
    
    if(iA==0){
      pvp.update_StateVariables(pc);  
      pvp.compute_right_vector(R, pc);
      norm_R = norm_R_0 = sqrt(R(i)*R(i));
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
  
    du(i) = -KI(i,j)*R(j);
  
    for(int ia=0; ia<DIM_3x3; ++ia)
      dM.data[ia] = du.get(ia);

    dpc = du.get(DIM_3x3);      

    Tensor<2> temp_M = pvp.sv.M(i,j);
    double temp_pc   = pc;  
/*
    if(fabs(dpc)>mat->cf_pcinf){
      Tensor<2> pFnp1, MI, M;
      M = pvp.sv.M(i,j) + dM(i,j);
      err += inv(M, MI);
      TensorA<2> pFn(pvp.Fs.pFn);            
      pFnp1(i,j) = MI(i,k)*pFn(k,l)*pvp.N(l,j);      
      double pJ = ttl::det(pFnp1);
      double temp_pc_np1 = pvp.compute_pc(pc, pJ);
      dpc = temp_pc_np1 - pc;
      du.get(DIM_3x3) = dpc;
      printf("pc_np1(%e) = pc_n(%e) + dpc(%e) \n", pc + dpc, pc, dpc);
    }
 */         
    double weight = 1.0;

    for(int iB=0; iB<64; iB++){
      pvp.sv.M(i,j) = temp_M(i,j) + weight*dM(i,j);
      pc = temp_pc + weight*dpc;
      
      pvp.update_StateVariables(pc);
      pvp.compute_right_vector(R, pc);

      double temp_norm_R = sqrt(R(i)*R(i));
      if(temp_norm_R<norm_R){        
        norm_R = temp_norm_R;
        break;
      } else {
        weight = weight/2.0;
        //printf("reduced weight (%e) is applied |R(try)| = %e, |R| = %e\n", weight, temp_norm_R, norm_R);
      }
    }  
  
 //   if(iA==pvp.solver_info->max_itr_M-1)
    //printf("(%d/%d)residual: |R| = %e, |R0| = %e, |R/R0| = %e\n", iA, pvp.solver_info->max_itr_M, norm_R, norm_R_0, norm_R/norm_R_0);

    if(norm_R/norm_R_0<pvp.solver_info->tol_M)
    {
      printf("(%d/%d)residual: |R| = %e, |R0| = %e, |R/R0| = %e\n", iA, pvp.solver_info->max_itr_M, norm_R, norm_R_0, norm_R/norm_R_0);  
      break;
    }
    
    double eng_norm = 0.0;
    for(int ia=0; ia<=DIM_3x3; ia++)
      eng_norm += weight*R.data[ia]*du.data[ia];
    
    eng_norm = fabs(eng_norm);

    // set initial energy norm
    if(iA==0)
      eng_norm_0 = eng_norm;

    if(eng_norm_0<pvp.solver_info->computer_zero)
      eng_norm_0 = pvp.solver_info->computer_zero;
      
    if(eng_norm/eng_norm_0 < (pvp.solver_info->tol_M)*(pvp.solver_info->tol_M))
    {
      if(1)
        printf("converge on energe norm %e\n", eng_norm/eng_norm_0);

     // break;
    }
    
  } 
  *pc_np1 = pc;
  
  return err; 
}

int poro_visco_plasticity_integration_algorithm_implicit(const MaterialPoroViscoPlasticity *mat,
                                                         const GcmSolverInfo *solver_info,
                                                         double *Fnp1,
                                                         double *Fn,
                                                         double *pFnp1,
                                                         double *pFn,
                                                         double *pc_np1,
                                                         const double pc_n)
{
  int err = 0;
  
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(mat);
  
  double dHdp = pvp.compute_dHdp(pc_n);
  if(fabs(dHdp)<1.0e-6)
  {
    //printf("dHdp = %.17e\n", dHdp);  
    return poro_visco_plasticity_integration_algorithm_staggered(mat, solver_info, Fnp1, Fn, pFnp1, pFn, pc_np1, pc_n);
  }
  pvp.set_tenosrs(Fnp1, Fn, pFnp1, pFn);
  pvp.set_scalars(*pc_np1, pc_n);
  pvp.set_solver_info(solver_info);
  
  Matrix_10x10<1> R, du;
  Matrix_10x10<2> K = {};
  
  Tensor<2> dM;
  
  double pc = pc_n;
  double dpc = 0.0;
  
  double norm_R_0 = pvp.solver_info->computer_zero;
  double norm_R   = pvp.solver_info->computer_zero;
  double eng_norm_0 = pvp.solver_info->computer_zero;
        
  for(int iA = 0; iA<pvp.solver_info->max_itr_M; iA++){
    
    if(iA==0){
      pvp.update_StateVariables(pc);  
      pvp.compute_right_vector(R, pc);
      norm_R = norm_R_0 = sqrt(R(i)*R(i));
    }
  
    pvp.compute_left_matrix(K, pc);
  
    try{
      du = ttl::solve(K,R);
    }
    catch(const int solve_err){  // no solution
      if(1)
        printf( "Matrix is singular. The solution (Poro-visco-plasticity) could not be computed.\n");

      ++err;
      break;
    } 
  
    for(int ia=0; ia<DIM_3x3; ++ia)
      dM.data[ia] = -du.get(ia);

    dpc = -du.get(DIM_3x3);      
    
    pvp.sv.M(i,j) += dM(i,j);
    pc += dpc;
    pvp.update_StateVariables(pc);
    pvp.compute_right_vector(R, pc);

    double temp_norm_R = sqrt(R(i)*R(i));
    norm_R = temp_norm_R;

    if(solver_info->debug)
      printf("(%d/%d)residual: |R| = %e, |R0| = %e, |R/R0| = %e\n", iA, pvp.solver_info->max_itr_M, norm_R, norm_R_0, norm_R/norm_R_0);

    if(norm_R/norm_R_0<pvp.solver_info->tol_M)
      break;
    
    double eng_norm = 0.0;
    for(int ia=0; ia<=DIM_3x3; ia++)
      eng_norm += R.data[ia]*du.data[ia];
    
    eng_norm = fabs(eng_norm);

    // set initial energy norm
    if(iA==0)
      eng_norm_0 = eng_norm;

    if(eng_norm_0<pvp.solver_info->computer_zero)
      eng_norm_0 = pvp.solver_info->computer_zero;
      
    if(eng_norm/eng_norm_0 < (pvp.solver_info->tol_M)*(pvp.solver_info->tol_M))
    {
      if(0)
        printf("converge on energe norm %e\n", eng_norm/eng_norm_0);

      //break;
    }
    
  } 
  *pc_np1 = pc;
  
  return err; 
}

int poro_visco_plasticity_integration_algorithm_implicit_yy(const MaterialPoroViscoPlasticity *mat,
                                                            const GcmSolverInfo *solver_info,
                                                            double *Fnp1,
                                                            double *Fn,
                                                            double *pFnp1,
                                                            double *pFn,
                                                            double *pc_np1,
                                                            const double pc_n)
{
  int err = 0;
  
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(mat);
  pvp.set_tenosrs(Fnp1, Fn, pFnp1, pFn);
  pvp.set_scalars(*pc_np1, pc_n);
  
  Matrix_10x10<1> R, du;
  Matrix_10x10<2> K = {}, KI = {};
  
  Tensor<2> dM;
  
  Tensor<2> LpFnp1, LFnI;
  TensorA<2> LpFn(pFn), LFn(Fn), LFnp1(Fnp1);
  
  LFnI(i,j) = ttl::inverse(LFn)(i,j);
  LpFnp1(i,j) = LpFn(i,l)*LFnI(l,k)*LFnp1(k,j); //eFnI = LpFn*LFnI

  double pJ = ttl::det(LpFnp1);
  double pc = pc_n;//pvp.compute_pc(pc_n, pJ);        
  double dpc = 0.0;
  
  double norm_R_0 = pvp.solver_info->computer_zero;
  double norm_R   = pvp.solver_info->computer_zero;
        
  for(int iA = 0; iA<pvp.solver_info->max_itr_M; iA++){
    
    if(iA==0){
      pvp.update_StateVariables(pc);  
      pvp.compute_right_vector(R, pc);
      norm_R = norm_R_0 = sqrt(R(i)*R(i));
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
  
    du(i) = -KI(i,j)*R(j);
  
    for(int ia=0; ia<DIM_3x3; ++ia)
      dM.data[ia] = du.get(ia);

    dpc = du.get(DIM_3x3);

    pvp.sv.M(i,j) += dM(i,j);
    pc += dpc;
      
    pvp.update_StateVariables(pc);
    pvp.compute_right_vector(R, pc);

    double norm_R = sqrt(R(i)*R(i));
  
    printf("(%d/%d)residual: |R| = %e, |R0| = %e, |R/R0| = %e\n", iA, pvp.solver_info->max_itr_M, norm_R, norm_R_0, norm_R/norm_R_0);

    if(norm_R/norm_R_0<pvp.solver_info->tol_M)
      break;
  } 
  *pc_np1 = pc;
  
  return err; 
}

int poro_visco_plasticity_integration_algorithm_explicit(const MaterialPoroViscoPlasticity *mat,
                                                         const GcmSolverInfo *solver_info,
                                                         double *Fnp1,
                                                         double *Fn,
                                                         double *pFnp1,
                                                         double *pFn,
                                                         double *pc_np1,
                                                         const double pc_n)
                                                         
{  
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(mat);
  pvp.set_tenosrs(Fnp1, Fn, pFnp1, pFn);
  pvp.set_scalars(*pc_np1, pc_n);
  pvp.set_solver_info(solver_info);
  
  pvp.update_StateVariables(pc_n);
  Tensor<2> pD;
  TensorA<2> pF(pFnp1), pF_n(pFn);
  
  pvp.compute_pD(pD, pc_n);
  pF(i,j) = pF_n(i,j) + solver_info->dt*pD(i,k)*pF_n(k,j);
  double pJ = ttl::det(pF);
  *pc_np1 = pvp.compute_pc(pc_n, pJ);
  
  return 0;
}

int poro_visco_plasticity_integration_algorithm(const MaterialPoroViscoPlasticity *mat,
                                                const GcmSolverInfo *solver_info,
                                                double *Fnp1,
                                                double *Fn,
                                                double *pFnp1,
                                                double *pFn,
                                                double *pc_np1,
                                                const double pc_n,
                                                const bool is_implicit)
{
  if(is_implicit)
    return poro_visco_plasticity_integration_algorithm_implicit(mat, solver_info, Fnp1, Fn, pFnp1, pFn, pc_np1, pc_n);
  else
    return poro_visco_plasticity_integration_algorithm_explicit(mat, solver_info, Fnp1, Fn, pFnp1, pFn, pc_np1, pc_n);
}

/// compute conforming pressure for given plastic deformation 
/// \param[in] pJ        determinant of pF
/// \param[in] mat_pvp   poro_viscoplaticity material object
/// \return    computed pc value  
double poro_visco_plasticity_compute_pc(double pJ, 
                                        double pc,
                                        const MaterialPoroViscoPlasticity *mat,
                                        const GcmSolverInfo *solver_info)
{
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(mat);
  pvp.set_solver_info(solver_info);
        
  return pvp.compute_pc(pc, pJ);
}



void poro_visco_plasticity_update_elasticity(double *eS_out,
                                             double *L_out,
                                             const MaterialPoroViscoPlasticity *param,
                                             double *eF_in,
                                             const double pc,
                                             const bool compute_4th_order){
  TensorA<2> eF(eF_in), eS(eS_out);
  TensorA<4> L(L_out);
  
  PvpMaterial mat(param);      
  double c  = mat.compute_c(pc);
  double d  = mat.compute_d(pc);
  double K  = mat.compute_bulk_modulus(c, d);
  double mu = mat.compute_shear_modulus(c, d);
  double coeff_U_alpha = mat.compute_coeff_U_alpha(c);
  double coeff_U_beta  = mat.compute_coeff_U_beta(c);

  PvpElasticity elasticity;
  elasticity.compute_elasticity_tensor(eS, L, eF, mu, K, c, coeff_U_alpha, 0.0, coeff_U_beta, compute_4th_order);
}

void poro_visco_plasticity_update_elasticity_dev(double *eS_out,
                                                 double *L_out,
                                                 const MaterialPoroViscoPlasticity *param,
                                                 double *eF_in,
                                                 const double pc,
                                                 const bool compute_4th_order){
  TensorA<2> eF(eF_in), eS(eS_out);
  TensorA<4> L(L_out);
  
  PvpMaterial mat(param);      
  double c  = mat.compute_c(pc);
  double d  = mat.compute_d(pc);
  double mu = mat.compute_shear_modulus(c, d);

  PvpElasticity elasticity;

  Tensor<2> dWdC, eC;
  eC(i,j) = eF(k,i)*eF(k,j);
  elasticity.compute_dWdC_dev(dWdC, L, eC, mu, compute_4th_order);

  eS(i,j) = 2.0*dWdC(i,j);
}

double poro_visco_plasticity_intf_compute_dudj(const double eJ,
                                               const double pc,
                                               const MaterialPoroViscoPlasticity *param){  
  PvpMaterial mat(param);      
  double c  = mat.compute_c(pc);
  double d  = mat.compute_d(pc);
  double K  = mat.compute_bulk_modulus(c, d);
  double coeff_U_alpha = mat.compute_coeff_U_alpha(c);
  double coeff_U_beta  = mat.compute_coeff_U_beta(c);

  PvpElasticity elasticity;
  return elasticity.compute_dUdJ(eJ, K, c, coeff_U_alpha, 0.0, coeff_U_beta);
}

double poro_visco_plasticity_intf_compute_d2udj2(const double eJ,
                                                 const double pc,
                                                 const MaterialPoroViscoPlasticity *param){  
  PvpMaterial mat(param);      
  double c  = mat.compute_c(pc);
  double d  = mat.compute_d(pc);
  double K  = mat.compute_bulk_modulus(c, d);
  double coeff_U_alpha = mat.compute_coeff_U_alpha(c);
  double coeff_U_beta  = mat.compute_coeff_U_beta(c);

  PvpElasticity elasticity;
  return elasticity.compute_d2UdJ2(eJ, K, c, coeff_U_alpha, 0.0, coeff_U_beta);
}

double poro_visco_plasticity_hardening(const double pc,
                                       const MaterialPoroViscoPlasticity *param){
  PvpMaterial mat(param);
  return mat.compute_H(pc);
}

void poro_visco_plasticity_intf_compute_gammas(double &gamma_dot_d,
                                               double &gamma_dot_v,
                                               double *Fnp1_in,
                                               double *Fn_in,
                                               double *pFnp1_in,
                                               double *pFn_in,
                                               const double pc_np1,
                                               const double pc_n,
                                               const double dt,
                                               const MaterialPoroViscoPlasticity *param,
                                               const GcmSolverInfo *solver_info)
{
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(param);       
  pvp.set_tenosrs(Fnp1_in, Fn_in, pFnp1_in, pFn_in); 

  pvp.set_scalars(pc_np1, pc_n);
  pvp.set_solver_info(solver_info);
  pvp.update_StateVariables(pc_np1);
    
  double d = pvp.mat.compute_d(pc_np1);
  gamma_dot_d = pvp.mat.compute_gamma_dot_d(d, pvp.sv.bar_tau, pvp.sv.g_tau);
  gamma_dot_v = pvp.mat.compute_gamma_dot_v(pvp.sv.pi, pvp.sv.pi_m, pvp.sv.g_pi);  
}


