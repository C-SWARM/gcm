/**
 * Authors:
 *  Sangmin Lee, [1], <slee43@nd.edu>
 *  [1] - University of Notre Dame, Notre Dame, IN
 */

#include"constitutive_model.h"
#include"poro_visco_plasticity.h"
#include"GcmSolverInfo.h"

constexpr const int Err = 1;
constexpr const double one_third = 1.0/3.0;

Tensor<2> I  = ttl::identity(i,j);
Tensor<2> O  = ttl::identity(i,j)*0.0;
      
template<int R, int D = 10, class S = double>
using Matrix_10x10  = ttl::Tensor<R, D, S>;

double Macaulay(const double x,
                const double x0){
  return (x >= x0)? x-x0 : 0.0;  
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

/// define material class for PVP model. This object is used only in this file scope
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
      return (param->p0 + c)*(d-1.0/d)/param->kappa;
    }

    /// compute material constant for volumetric part of strain energy density function
    /// U = 1/2*K*(J-1)*ln(J) + c*ln(J) + (alpha1 + alpha2*ln(J))*exp(-beta*ln(J))
    /// note: alpha2 = 0 for computing stress and elasticity tensor
    ///
    /// \param[in] c pvp model function value of pc
    /// \return    material constant (alpha_1)
    double compute_coeff_U_alpha(const double c){
      return (param->p0 + c)*pow(param->c_inf/(param->c_inf - c), param->n)*param->kappa;
    }

    /// compute material constant for volumetric part of strain energy density function
    /// U = 1/2*K*(J-1)*ln(J) + c*ln(J) + (alpha1 + alpha2*ln(J))*exp(-beta*ln(J))
    /// note: alpha2 = 0 for computing stress and elasticity tensor
    ///
    /// \param[in] c pvp model function value of pc
    /// \return    material constant (beta)        
    double compute_coeff_U_beta(const double c){
      return pow((param->c_inf - c)/param->c_inf, param->n)/param->kappa;
    } 
    
    /// compute a
    ///
    /// \param[in] pc
    /// \param[in] c pvp model function value of pc
    /// \return    a
    double compute_a(const double pc,
                     const double c){
      return (pc + c)/(1.0 + param->alpha);
    }
    
    /// compute c
    ///
    /// \param[in] pc
    /// \return    c
    double compute_c(const double pc){
      return param->c_inf*(1.0-exp(-param->Gamma*Macaulay(pc, param->pc_b)));
    }

    /// compute b
    ///
    /// \param[in] pi   Kirchhoff pressure(-tr(tau)/3.0)
    /// \param[in] pi_m parameter for yeild surface 
    /// \param[in] c    pvp model function value of pc
    /// \return    b   
    double compute_b(const double pi,
                     const double pi_m){
      return (pi >= pi_m)? param->alpha : 1.0;
    }

    /// compute d
    ///
    /// \param[in] pc
    /// \return    d    
    double compute_d(const double pc){
      return 1.0 + param->B*Macaulay(pc, param->pc_b);
    }

    /// compute g_tau
    ///
    /// \param[in] pc
    /// \param[in] a pvp model function value of pc
    /// \return    g_tau          
    double compute_g_tau(const double pc,
                         const double a){
      return sqrt(1.5)*a*param->M; // sqrt(3.0/2.0) = sqrt(1.5)
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
      return -param->a1*exp(-param->Lambda1/pc) - param->a2*exp(-param->Lambda2/pc);
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
      return param->gamma_dot_0*(1.0-1.0/d)*pow(bar_tau/g_tau, 1.0/param->m);
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
      return (pi>pi_m)? param->gamma_dot_0*pow((pi-pi_m)/g_pi, 1.0/param->m): 0.0;
    }

    /// compute dilatancy function value. zero for this model 
    /// \return 0.0 always 
    double compute_beta_D(void){return 0.0;}

    /// compute compaction function value 
    /// 
    /// \param[in] pc
    /// \return    beta_C     
    double compute_beta_C(const double pc){
      return param->g0*(1.0 - pc/param->pc_inf);
    }
    
    /// compute derivative of H w.r.t p_c
    ///
    /// \param[in] pc
    /// \return    dHdp
    double compute_dHdp(const double pc){
      return -1.0/(pc*pc)*(param->a1*param->Lambda1*exp(-param->Lambda1/pc) + param->a2*param->Lambda2*exp(-param->Lambda2/pc));
    }
    
    /// compute derivative of c(p_c) w.r.t p_c
    ///
    /// \param[in] pc
    /// \return    dcdp
    double compute_dcdp(const double pc){
      return param->c_inf*param->Gamma*exp(-param->Gamma*Macaulay(pc, param->pc_b));
    }

    /// compute derivative of d(p_c) w.r.t p_c
    ///
    /// \param[in] pc
    /// \return    dddp    
    double compute_dddp(const double pc){
      return (pc>param->pc_b)? param->B: 0.0;
    }

    /// compute derivative of a(p_c) w.r.t p_c
    ///
    /// \param[in] pc
    /// \return    dadp    
    double compute_dadp(const double pc){
      return 1.0/(1.0+param->alpha)*(1.0 + compute_dcdp(pc));
    }

    /// compute derivative of g_tau(p_c) w.r.t p_c
    ///
    /// \param[in] pc
    /// \return    dg_tau_dp    
    double compute_dg_tau_dp(const double pc){
      return sqrt(1.5)*compute_dadp(pc)*param->M; // sqrt(3.0/2.0) = sqrt(1.5)
    }            
};

/// This class carries all variables that need to be updated when 
/// pFnp1 or M is updated such as eS, eF, and tau are updated.
/// This is a local class.
class StateVariables
{
  public:
    Tensor<2> M, MI, eFnp1, tau, eS, eSd, hat_tau, psi_d, psi_v;;
    Tensor<4> L;
    double pJ, pi, bar_tau, g_tau, g_pi, pi_m;
    
    StateVariables(){};
};

/// object for the elasticity part of the pvp model.
/// This is a local class.
class PvpElasticity
{
  public:

    /// compute derivative of deviatoric part of W(strain energy density function) w.r.t eC
    ///
    /// \param[out] dWdC   computed derivative of W (deviatoric part, 2nd order tensor)
    /// \param[out] d2WdC2 computed 2nd derivative of W (deviatoric part, 4th order tensor)
    /// \param[in]  C      C = eF'eF (right Cauchy Green tensor)
    /// \param[in]  mu     shear modulus
    /// \param[in]  compute_4th_order if yes: compute d2WdC2, default is false
    ///                                  no : skip computing d2WdC2
    template<class T1, class T2, class T3> void compute_dWdC_dev(T1 &dWdC,
                                                                 T2 &d2WdC2,
                                                                 T3 &C,
                                                                 const double mu,
                                                                 bool compute_4th_order = false){

      double CJ = ttl::det(C); // eCJ = eJ*eJ; pow(eCJ, -1.0/3.0) = pow(eJ, -2.0/3.0)
      double factor = pow(CJ, -one_third);
      
      double trC = C(i,i);
      
      Tensor<2> CI;
      int err = inv(C, CI);
      
      dWdC(i,j) = 0.5*mu*factor*(I(i,j) - one_third*trC*CI(i,j));
    
      if(compute_4th_order)
      {
        Tensor<4> CIxCI, dCIdC;
        CIxCI(i,j,k,l) = CI(i,j)*CI(k,l);
        dCIdC(i,j,k,l)  = -CI(i,k)*CI(l,j);
        
        d2WdC2(i,j,k,l) = 0.5*mu*factor*(-one_third*CI(i,j)*I(k,l) - one_third*I(i,j)*CI(k,l) 
                                         + one_third*one_third*trC*CIxCI(i,j,k,l) - one_third*trC*dCIdC(i,j,k,l));
      }
      if(err>0)
        throw Err;
        
    }
    
    /// compute derivative of volumetric part of W(strain energy density function, U) w.r.t eJ
    /// U = 1/2*K*(J-1)*ln(J) + c*ln(J) + (alpha1 + alpha2*ln(J))*exp(-beta*ln(J))
    /// 
    /// \param[in] J      det(eF)
    /// \param[in] K      bulk modulus
    /// \param[in] c      material parameter for pvp as a functin of p_c
    /// \param[in] alpha1 p_c dependent contribution to particle compaction
    /// \param[in] alpha2 p_c dependent contribution to particle compaction
    /// \param[in] beta   p_c dependent contribution to particle compaction
    /// \return    computed dUdJ
    double compute_dUdJ(const double J,
                        const double K,
                        const double c,
                        const double alpha1,
                        const double alpha2,
                        const double beta){
      return 0.5*K*(log(J) + (J-1.0)/J) + c/J + (alpha2 - beta*(alpha1 + alpha2*log(J)))*exp(-beta*log(J))/J;
    }

    /// compute 2nd derivative of volumetric part of W(strain energy density function, U) w.r.t eJ
    /// U = 1/2*K*(J-1)*ln(J) + c*ln(J) + (alpha1 + alpha2*ln(J))*exp(-beta*ln(J))
    /// 
    /// \param[in] J      det(eF)
    /// \param[in] K      bulk modulus
    /// \param[in] c      material parameter for pvp as a functin of p_c
    /// \param[in] alpha1 p_c dependent contribution to particle compaction
    /// \param[in] alpha2 p_c dependent contribution to particle compaction
    /// \param[in] beta   p_c dependent contribution to particle compaction
    /// \return    computed d2UdJ2    
    double compute_d2UdJ2(const double J,
                          const double K,
                          const double c,
                          const double alpha1,
                          const double alpha2,
                          const double beta){
      double JJ = 1.0/(J*J);
      return 0.5*K*(1.0/J + 1.0*JJ) - c*JJ + ((1.0+beta)*beta*(alpha1 + alpha2*log(J))
                                                 -(1.0+2.0*beta)*alpha2)*exp(-beta*log(J))*JJ;
    }

    /// compute PKII and elasticity tensor w.r.t eC
    /// 
    /// \param[out] eS     PKII (2nd order tensor)
    /// \param[out] L      elasticity tensor (4th order tensor)   
    /// \param[in]  eF     elastic part of deformation gradient tensor
    /// \param[in]  mu     shear modulus
    /// \param[in]  K      bulk modulus
    /// \param[in]  c      material parameter for pvp as a functin of p_c
    /// \param[in]  alpha1 p_c dependent contribution to particle compaction
    /// \param[in]  alpha2 p_c dependent contribution to particle compaction
    /// \param[in]  beta   p_c dependent contribution to particle compaction
    /// \param[in]  compute_4th_order if yes: compute L, default is false
    ///                                  no : skip computing L    
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
        L(i,j,k,l) = 4.0*L(i,j,k,l) + (eJ*dUdJ + eJ*eJ*d2UdJ2)*eCI(i,j)*eCI(k,l) - 2.0*eJ*dUdJ*eCI(i,k)*eCI(l,j);
      }
      
      if(err>0)
        throw Err;
    }
};

/// All computation for PVP model is done through this class.
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
    double dt;  
            
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
    /// \param[in] total_Lagrangian if true, compute tensors based on total lagrangian, default is false (updated Lagrangian)
    void set_tenosrs(double *Fnp1,
                     double *Fn,
                     double *pFnp1,
                     double *pFn,
                     double *hFnp1,
                     double *hFn,
                     bool total_Lagrangian = false){
      Fs.Fnp1  = Fnp1;
      Fs.pFnp1 = pFnp1;
      Fs.hFnp1 = hFnp1;
 
      if(total_Lagrangian){
        Fs.Fn    = I.data;
        Fs.pFn   = I.data;
        Fs.hFn   = I.data;
      } else {
        Fs.Fn    = Fn;
        Fs.pFn   = pFn;
        Fs.hFn   = hFn;
      }
      compute_tensors();
    }
    
    /// senumerical parameters such as tolerance and maximum number of iterations
    ///
    /// \param[in] max_itr_in       maximum number of NR iterations
    /// \param[in] tol_in           NR tolerance
    /// \param[in] computer_zero_in computer zero
    void set_solver_info(const GcmSolverInfo *solver_info_in,
                         const double dt_in){
      solver_info = solver_info_in;
      dt          = dt_in;
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
    /// \param[in] update4integration if true : update M first and compute other dependent variables
    ///                                  false: all variables are computed based on given pFnp1 from PDE
    void update_StateVariables(const double pc,
                               bool update4integration = true){
      
      int err = 0;
      
      TensorA<2> pFnp1(Fs.pFnp1);
      if(update4integration){

        err += inv(sv.M, sv.MI);
        TensorA<2> pFn(Fs.pFn), pFnp1(Fs.pFnp1);
        pFnp1(i,j) = sv.MI(i,k)*pFn(k,l)*N(l,j);

        sv.eFnp1(i,j) = Fa(i,k)*sv.M(k,j);
      
      } else {

        Tensor<2> FaI, pFnp1I, hFnp1I;
        TensorA<2> Fnp1(Fs.Fnp1), hFnp1(Fs.hFnp1);
  
        err += inv(pFnp1, pFnp1I);
        err += inv(Fa,    FaI);
        err += inv(hFnp1, hFnp1I);

        sv.eFnp1(i,j) = Fnp1(i,k)*hFnp1I(k,l)*pFnp1I(l,j);  
        sv.M(i,j) = FaI(i,k)*sv.eFnp1(k,j); // pvp.set_tenosrs dosen't update M correcltly (M(i,j) = eFnI(i,k)*FrI(k,l)*eFn(l,j);)
        err += inv(sv.M, sv.MI);
      }
      sv.pJ = ttl::det(pFnp1);
     
      double c  = mat.compute_c(pc);
      double d  = mat.compute_d(pc);
      double K  = mat.compute_bulk_modulus(c, d);
      double mu = mat.compute_shear_modulus(c, d);

      double coeff_U_alpha = mat.compute_coeff_U_alpha(c);
      double coeff_U_beta  = mat.compute_coeff_U_beta(c);

      elasticity.compute_elasticity_tensor(sv.eS, sv.L, sv.eFnp1, mu, K, c, coeff_U_alpha, 0.0, coeff_U_beta, true);
      
      double tr_eS = sv.eS(i,i);
      sv.eSd(i,j) = sv.eS(i,j) - one_third*tr_eS*I(i,j);
      sv.tau(i,j) = sv.pJ*sv.eFnp1(i,k)*sv.eS(k,l)*sv.eFnp1(j,l);
      sv.pi       = compute_pi(sv.tau);
      
      sv.hat_tau(i,j)  = sv.tau(i,j) + sv.pi*I(i,j);
      sv.bar_tau = sqrt(sv.hat_tau(i,j)*sv.hat_tau(i,j));
      
      double a = mat.compute_a(pc, c);
      sv.pi_m  = mat.compute_pi_m(a, c);
      double b = mat.compute_b(sv.pi, sv.pi_m);
      
      sv.g_tau = mat.compute_g_tau(pc, a);
      sv.g_pi  = mat.compute_g_pi(a, b);
      
      compute_psi_d(sv.psi_d);
      compute_psi_v(sv.psi_v, pc);            
    }
    
    /// compute Krichhoff pressure, pi = tr(tau)
    /// internal vaiables define in sv(StateVariables) class are used.
    ///
    /// \param[in] tau tau = pJFeSF'
    /// \return    computed pi
    template<class T> double compute_pi(const T &tau){
      return -tau(i,i)*one_third;
    }

    /// compute psi_d
    ///
    /// \param[out] psi_d = hat(eS)/||hat(eS)||
    /// \return     computed pi    
    template<class T> void compute_psi_d(T &psi_d){

      double m_eSd = sqrt(sv.eSd(i,j)*sv.eSd(i,j));
      double beta_D = mat.compute_beta_D();
      
      if(m_eSd<solver_info->computer_zero)
        psi_d(i,j) = one_third*beta_D*I(i,j);
      else
        psi_d(i,j) = sv.eSd(i,j)/m_eSd + one_third*beta_D*I(i,j);
    }

    /// compute psi_v
    ///
    /// \param[out] psi_v = 1/3 beta_C delta_ij
    /// \param[in]  pc conforming pressure
    template<class T> void compute_psi_v(T &psi_v,
                                         const double pc){
      double beta_C = mat.compute_beta_C(pc);
      psi_v(i,j) = -one_third*beta_C*I(i,j);
    }
    
    /// compute_residual RM
    ///
    /// \param[out] RM computed residual for M
    /// \param[in]  pc conforming pressure
    template<class T> void compute_RM(T &RM,
                                      const double pc){
      Tensor<2> pD;
      compute_pD(pD, pc);      
      RM(i,j) = dt*pD(i,j) - I(i,j) + A(i,k)*sv.M(k,j);
    }
    
    /// compute residual Rpc
    /// 
    /// \param[out] Rpc computed residual for p_c
    /// \param[in]  pc conforming pressure    
    void compute_Rp(double &Rpc,
                    const double pc){
                      
      Rpc = sv.pJ - exp(mat.compute_H(pc));
    }
    
    /// compute d_bar_tau_d_tau (d_bar_tau_d_tau = d_bar_tau/d||tau||
    ///
    /// \param[out] d_bar_tau computed d_bar_tau/d||tau||
    template<class T> void compute_d_bar_tau_d_tau(T &d_bar_tau){

      if(sv.bar_tau<solver_info->computer_zero)
        d_bar_tau = {};
      else
        d_bar_tau(i,j) = sv.hat_tau(i,j)/sv.bar_tau;      
    }
    
    /// compute dS/dM
    ///
    /// \param[out] dSdM computed dS/dM
    template<class T> void compute_dSdM(T &dSdM){
      Tensor<2> FaTFa;
      
      FaTFa(i,j) = Fa(k,i)*Fa(k,j);                                  
      dSdM(i,j,k,l) = 0.5*sv.L(i,j,m,r)*(I(m,l)*FaTFa(k,v)*sv.M(v,r) + sv.M(v,m)*FaTFa(v,k)*I(r,l));
    }
    
    /// compute dtau/dM
    ///
    /// \param[out] d_tau_dM computed dtau/dM
    /// \param[in]  dSdM     dS/dM computed prior to this function
    ///                      because dS/dM is used multiple places for the pvp integration.
    template<class T1, class T2> void compute_d_tau_dM(T1 &d_tau_dM,
                                                       const T2 &dSdM){
      Tensor<4> deF_dM, deFT_dM;     

      deF_dM(i,j,k,l)  = Fa(i,k)*I(j,l);
                  
      d_tau_dM(i,j,k,l) = -sv.tau(i,j)*sv.MI(l,k) + sv.pJ*deF_dM(i,m,k,l)*sv.eS(m,r)*sv.eFnp1(j,r)
                                                  + sv.pJ*sv.eFnp1(i,m)*dSdM(m,r,k,l)*sv.eFnp1(j,r)
                                                  + sv.pJ*sv.eFnp1(i,m)*sv.eS(m,r)*deF_dM(j,r,k,l);
    }
    
    /// compute dgamma_dot_d/dM
    /// 
    /// \param[out] dgamma_dot_d_dM computed dgamma_dot_d/dM
    /// \param[in]  d_tau_dM        computed prior since it is used multiple places
    /// \param[in]  d               d(p_c)
    template<class T> void compute_dgamma_dot_d_dM(T &dgamma_dot_d_dM,
                                                   const Tensor<4> &d_tau_dM,
                                                   const double d){
      double factor = mat.param->gamma_dot_0/(mat.param->m*sv.g_tau)*(1.0 - 1.0/d)*pow(sv.bar_tau/sv.g_tau, 1.0/mat.param->m-1.0);
      Tensor<2> d_bar_tau;
      compute_d_bar_tau_d_tau(d_bar_tau);
      
      dgamma_dot_d_dM(k,l) = factor*d_bar_tau(i,j)*d_tau_dM(i,j,k,l);
    }
    
    /// compute dgamma_dot_v/dM
    ///
    /// \param[out] dgamma_dot_v/dM computed dgamma_dot_v/dM
    /// \param[in]  d_tau_dM        computed prior since it is used multiple places
    template<class T> void compute_dgamma_dot_v_dM(T &dgamma_dot_v_dM,
                                                   const Tensor<4> &d_tau_dM){
      if(sv.pi>sv.pi_m){        
        double factor = -mat.param->gamma_dot_0/(3.0*mat.param->m*sv.g_pi)*pow((sv.pi-sv.pi_m)/sv.g_pi, 1.0/mat.param->m - 1.0);
        dgamma_dot_v_dM(k,l) = factor*I(i,j)*d_tau_dM(i,j,k,l);
      }
      else
        dgamma_dot_v_dM(k,l) = O(k,l);

    }
    
    /// compute dpsi_d/dM
    /// 
    /// \param[in] dpsi_d_dM computed dpsi_d/dM
    /// \param[in]  dSdM     dS/dM computed prior to this function, it is used multiple places.
    template<class T1, class T2> void compute_dpsi_d_dM(T1 &dpsi_d_dM,
                                                        const T2 &dSdM){
      double m_eSd = sqrt(sv.eSd(i,j)*sv.eSd(i,j));
      
      if(m_eSd<solver_info->computer_zero)
        dpsi_d_dM = {};
      else
        dpsi_d_dM(i,j,k,l) = 1.0/m_eSd*(dSdM(i,j,k,l) - I(i,j)*I(m, r)*dSdM(m,r,k,l)*one_third)
                            -1.0/m_eSd/m_eSd/m_eSd*sv.eSd(i,j)*sv.eSd(m,r)*dSdM(m,r,k,l);
    }
    
    /// compute tangent for RM w.r.t M
    ///
    /// \param[out] dRMdM computed tangent for RM w.r.t M
    /// \param[in]  dSdM  dS/dM computed prior to this function, it is used multiple places.
    /// \param[in]  pc    conforming pressure    
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
      
      dRMdM(i,j,k,l) = A(i,k)*I(j,l) + dt*(sv.psi_d(i,j)*dgamma_dot_d_dM(k,l) + gamma_dot_d*dpsi_d_dM(i,j,k,l) 
                                            + sv.psi_v(i,j)*dgamma_dot_v_dM(k,l));
    }
   
    /// compute p_c dependent PKII from dW/dp w.r.t eC
    ///
    /// \param[in] dE_dWdp computed PKII
    /// \param[in]  pc    conforming pressure 
    template<class T> void compute_dE_dWdp(T &dE_dWdp,
                                           const double pc){
      double c = mat.compute_c(pc);
      double d = mat.compute_d(pc);
            
      double dcdp = mat.compute_dcdp(pc);
      double ddpc = mat.compute_dddp(pc);
      double dmudp = (dcdp*(d-1.0/d) + c*(1.0 + 1.0/(d*d))*ddpc)*mat.param->mu_1;
      double dKdp  = (dcdp*(d-1.0/d) + (mat.param->p0 + c)*(1.0 + 1.0/(d*d))*ddpc)/mat.param->kappa;
      
      double pow_of_c = pow(mat.param->c_inf/(mat.param->c_inf - c), mat.param->n);

      double coeff_U_alpha1 = dcdp*pow_of_c*mat.param->kappa*(1.0 + mat.param->n*(mat.param->p0 + c)/(mat.param->c_inf - c));
      double coeff_U_alpha2 = dcdp*(mat.param->p0 + c)/(mat.param->c_inf - c)*mat.param->n;
            
      double coeff_U_beta  = mat.compute_coeff_U_beta(c);

      Tensor<4> L;
      elasticity.compute_elasticity_tensor(dE_dWdp, L , sv.eFnp1, dmudp, dKdp, dcdp, coeff_U_alpha1, coeff_U_alpha2, coeff_U_beta, false);
    }
    
    /// compute d_tau/dp
    ///
    /// \param[in]  dSdM  dS/dM computed prior to this function, it is used multiple places.
    /// \param[in]  pc    conforming pressure 
    template<class T1, class T2> void compute_d_tau_dp(T1 &d_tau_dp,
                                                       const T2 &dSdp,
                                                       const double pc){
      double dHdp = mat.compute_dHdp(pc);                                 
      d_tau_dp(i,j) = sv.eFnp1(i,k)*sv.pJ*dSdp(k,l)*sv.eFnp1(j,l) - dHdp*sv.tau(i,j);
    }
    
    /// compute d_bar_tau*d_tau/dp
    ///
    /// \param[out] d_tau_dp computed d_bar_tau*d_tau/dp
    template<class T> double compute_d_bar_tau_dp(T &d_tau_dp){
      Tensor<2> d_bar_tau;
      
      compute_d_bar_tau_d_tau(d_bar_tau);
      return d_bar_tau(i,j)*d_tau_dp(i,j);      
    }
    
    /// compute d_bar_tau/dp - d_tau/dp
    ///
    /// \param[in] d_tau_dp computed d_bar_tau/dp - d_tau/dp
    /// \param[in] pc       conforming pressure 
    template<class T> double compute_d_bar_tau_g_tau_dp(const T &d_tau_dp,
                                                        const double pc){
      double d_bar_tau_dp = compute_d_bar_tau_dp(d_tau_dp);
      
      double d_g_tau_dp = mat.compute_dg_tau_dp(pc);      
      return d_bar_tau_dp - sv.bar_tau/sv.g_tau*d_g_tau_dp;
    }

    /// compute tangent for RM w.r.t p_c
    ///
    /// \param[out] dRMdp computed tangent for RM w.r.t p_c
    /// \param[in]  pc    conforming pressure     
    template<class T> void compute_dRMdp(T &dRMdp,
                                         const double pc){
      // compute dSdp d_tau_dp
      double d = mat.compute_d(pc);
      
      Tensor<2> dSdp, d_tau_dp;
      compute_dE_dWdp(dSdp, pc);
      
      
      double dHdp = mat.compute_dHdp(pc);
      Tensor<2> eC;
      eC(i,j) = sv.eFnp1(k,i)*sv.eFnp1(k,j);      
      dSdp(i,j) = dSdp(i,j) - dHdp*sv.L(i,j,k,l)*eC(k,l); // effect from eFnp1 
            
      compute_d_tau_dp(d_tau_dp, dSdp, pc);
      
      double bar_tau_g_tau = sv.bar_tau/sv.g_tau;
      double pow_bar_tau_g_tau = pow(sv.bar_tau/sv.g_tau, 1.0/mat.param->m);      
                                                
      // compute dgamma_dot_d_pc
      double factor1 = mat.param->gamma_dot_0/(d*d)*pow_bar_tau_g_tau;
      
      double dddp = mat.compute_dddp(pc);      
      double d_bar_tau_g_tau_dp = compute_d_bar_tau_g_tau_dp(d_tau_dp,pc);

      double dgamma_dot_d_pc = factor1*dddp;      
      if(sv.bar_tau>solver_info->computer_zero)
        dgamma_dot_d_pc += mat.param->gamma_dot_0/(mat.param->m*sv.g_tau)*(1.0-1.0/d)*pow_bar_tau_g_tau/bar_tau_g_tau*d_bar_tau_g_tau_dp;
      
      // compute gamma_dot_d
      double gamma_dot_d = mat.compute_gamma_dot_d(d, sv.bar_tau, sv.g_tau);
      
      // compute d_psi_d_dp
      Tensor<2> d_psi_d_dp = {};
      double m_eSd = sqrt(sv.eSd(i,j)*sv.eSd(i,j));

      if(m_eSd>solver_info->computer_zero){        
        double tr_dSdp = dSdp(i,i);      
        double eSd_dSdp = 1.0/(m_eSd*m_eSd*m_eSd)*sv.eSd(i,j)*dSdp(i,j);
        d_psi_d_dp(i,j) = 1.0/m_eSd*(dSdp(i,j) - tr_dSdp*I(i,j)*one_third) - eSd_dSdp*sv.eSd(i,j);
      }
      
      // compute dgamma_dot_v_pc
      double factor3 = 0.0;
      if(sv.pi>sv.pi_m)
        factor3 = mat.param->gamma_dot_0/(mat.param->m*sv.g_pi)*pow((sv.pi-sv.pi_m)/sv.g_pi, 1.0/mat.param->m-1.0);

      double dadpc = mat.compute_dadp(pc);
      double dcdp = mat.compute_dcdp(pc);
      double dpidp = -d_tau_dp(i,i)*one_third;
      double dpi_mdp = dadpc - dcdp;
      double dg_pidp = mat.compute_b(sv.pi, sv.pi_m)*dadpc;
      double dgamma_dot_v_pc = factor3*(dpidp - dpi_mdp - (sv.pi-sv.pi_m)/sv.g_pi*dg_pidp);

      // compute gamma_dot_v
      double gamma_dot_v = mat.compute_gamma_dot_v(sv.pi, sv.pi_m, sv.g_pi);
      
      // compute d_psi_v_dp
      Tensor<2> d_psi_v_dp = mat.param->g0/(3.0*mat.param->pc_inf)*I(i,j);

      dRMdp(i,j) = dt*(dgamma_dot_d_pc*sv.psi_d(i,j) + gamma_dot_d*d_psi_d_dp(i,j) + dgamma_dot_v_pc*sv.psi_v(i,j) + gamma_dot_v*d_psi_v_dp(i,j));      
    }
    
    /// compute pD = gamma_dot_d*psi_d + gamma_dot_v*psi_v
    ///
    /// \param[out] pD computed gamma_dot_d*psi_d + gamma_dot_v*psi_v
    /// \param[in]  pc conforming pressure   
    template<class T> void compute_pD(T &pD,
                                      const double pc){
      double d = mat.compute_d(pc);
      double gamma_dot_d = mat.compute_gamma_dot_d(d, sv.bar_tau, sv.g_tau);
      double gamma_dot_v = mat.compute_gamma_dot_v(sv.pi, sv.pi_m, sv.g_pi);
      pD(i,j) = gamma_dot_d*sv.psi_d(i,j) + gamma_dot_v*sv.psi_v(i,j);                                          
    }
    
    /// compute tagents of RM and Rpc w.r.t M and p_c
    /// K = [dRM/dM, dRpc/dM; dRpc/dM, dRpc/dpc]
    ///
    /// \param[out] K  computed tangents
    /// \param[in]  pc conforming pressure
    void compute_left_matrix(Matrix_10x10<2> &K,
                             const double pc){
      Tensor<4> dRMdM, dSdM;
      Tensor<2> dRMdp, dRpdM;
                  
      compute_dSdM(dSdM);
      compute_dRMdM(dRMdM, dSdM, pc);                    
      compute_dRMdp(dRMdp, pc);

      dRpdM(i,j) = -sv.pJ*sv.MI(j,i);    

      double dRpdp = -sv.pJ*mat.compute_dHdp(pc);
      
      for(int ia=0; ia<DIM_3x3; ia++){
        K(ia,DIM_3x3) = dRMdp.get(ia);
        K(DIM_3x3,ia) = dRpdM.get(ia);
        for(int ib=0; ib<DIM_3x3; ib++){
          K(ia,ib) = dRMdM.get(ia*DIM_3x3+ib);
        }
      }
      K(DIM_3x3,DIM_3x3) = dRpdp;
    }
    
    /// compute residuals, RM and Rpc
    /// R = [RM; dRpc]
    ///
    /// \param[out] R  computed residualss
    /// \param[in]  pc conforming pressure    
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
    
    /// compute pc as a function of pJ numerically
    ///
    /// \param[in] pcn conforming pressure at t(n)
    /// \param[in] pJ  det(pFnp1) 
    /// \return computed conforming pressure
    double compute_pc(double pcn,
                      const double pJ){
                        
      bool is_convg = true;
      
      double logJp = log(pJ);      
      double DfDpcn = mat.compute_dHdp(pcn);
      
      int it=0;
      double pcnp1 = pcn;
      
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
        if(it==maxIt)
          is_convg = false;
          
      } else {
        const int maxIt=100;
        pcnp1 = pcn;
        
        while(it < maxIt){
          ++it;
          double f = -logJp + mat.compute_H(pcnp1);
          double dHdp = mat.compute_dHdp(pcnp1);
          
          double dpc = -f/dHdp;
          pcnp1 = pcnp1 + dpc;          
          if(fabs(logJp - mat.compute_H(pcnp1)) < solver_info->tol_hardening)
            break;
        }
        if(it==maxIt)
          is_convg = false;        
      }
      
      if(!is_convg && solver_info->debug)
        printf("Cannot compute pc for pJ = %e: pc_n = %e, pc_np1 = %e\n", pJ, pcn, pcnp1);
        
      return pcnp1;                        
    }

    /// compute dMdU for an finite element for PDE
    /// 
    /// \param[out] dMdUs   computed dMdU, nne*ndofn number of 2nd order tensors
    /// \param[in]  pc      conforming pressure
    /// \param[in]  Grad_us Grad(delta u), nne*ndofn number of 2nd order tensors
    /// \param[in]  nne     number of nodes in an element
    /// \param[in]  ndofn   number of nodal degree of freedoms
    void compute_dMdu(double* dMdUs,
                      const double pc,
                      double* Grad_us,
                      const int nne,
                      const int ndofn){
                                          
      Tensor<4> U;
      Tensor<2> B, dE_dWdp = {}, eC;
      Tensor<4> II, chi, IoxI;

      II(i,j,k,l) = I(i,k)*I(j,l);
      IoxI(i,j,k,l) = I(i,j)*I(k,l);
      
      eC(i,j) = sv.eFnp1(k,i)*sv.eFnp1(k,j);
      
      double d = mat.compute_d(pc);
      double dddp = mat.compute_dddp(pc);
      double dHdp = mat.compute_dHdp(pc);
      double pi_pi_m_over_g_pi = 0.0;
      double pow_pi_pi_m_over_g_pi = 0.0;
      if(sv.pi>sv.pi_m)
      {  
        pi_pi_m_over_g_pi = (sv.pi-sv.pi_m)/sv.g_pi;
        pow_pi_pi_m_over_g_pi = pow(pi_pi_m_over_g_pi, 1.0/mat.param->m - 1.0);
      }
      
      // start computing A(left side);
      // compute chi
      Tensor<2> eFTFa;
      eFTFa(i,j) = sv.eFnp1(k,i)*Fa(k,j);
      
      compute_dE_dWdp(dE_dWdp, pc);
      chi(i,j,k,l) = 2.0*dHdp*sv.L(i,j,l,m)*eFTFa(m,k) - dE_dWdp(i,j)*sv.MI(l,k);
            
      //compute D(gamma_dot_d) part      
      double gamma_dot_d = mat.compute_gamma_dot_d(d, sv.bar_tau, sv.g_tau);
      double gamma_dot_v = mat.compute_gamma_dot_v(sv.pi, sv.pi_m, sv.g_pi);      

      double pow_bar_tau_g_tau = pow(sv.bar_tau/sv.g_tau, 1.0/mat.param->m);
      double d_g_tau_dp = mat.compute_dg_tau_dp(pc); 
      
      double factor1 = pow_bar_tau_g_tau*mat.param->gamma_dot_0*(-dddp/(d*d) 
                       + (1.0 - 1.0/d)/(mat.param->m*sv.g_tau)*d_g_tau_dp);
      double factor2 = sv.pJ*mat.param->gamma_dot_0/(mat.param->m*sv.g_tau)*(1.0 - 1.0/d)*pow_bar_tau_g_tau/sv.bar_tau*sv.g_tau;

      Tensor<2> U1, U2, d_bar_tau_d_tau, eFeSeFT, d_bar_tau_d_tau_sym;
      Tensor<2> eFTdtaueF;
      eFTdtaueF(i,j) = sv.eFnp1(k,i)*d_bar_tau_d_tau(k,l)*sv.eFnp1(l,j);
            
      compute_d_bar_tau_d_tau(d_bar_tau_d_tau);
      d_bar_tau_d_tau_sym(i,j) = 0.5*(d_bar_tau_d_tau(i,j) + d_bar_tau_d_tau(j,i));
      
      eFeSeFT(i,j) = sv.eFnp1(i,k)*sv.eS(k,l)*sv.eFnp1(j,l);
      
      double sub_factor1 = d_bar_tau_d_tau(i,j)*eFeSeFT(i,j);

      U1(i,j) = factor1*sv.MI(j,i);
      U2(i,j) = factor2*(-dHdp*sub_factor1*sv.MI(j,i) + dHdp*2.0*Fa(k,i)*d_bar_tau_d_tau_sym(k,l)*sv.eFnp1(l,m)*sv.eS(m,j)
                         + eFTdtaueF(k,l)*chi(k,l,i,j));
                         
      // compute D(psi_d) part
      Tensor<4> U3;      
      double m_eSd = sqrt(sv.eSd(i,j)*sv.eSd(i,j));

      Tensor<4> eSdoxeSd = {};
      
      if(m_eSd>solver_info->computer_zero){
        double one_over_m_eSd = 1.0/m_eSd;        
        eSdoxeSd(i,j,k,l) = one_over_m_eSd*one_over_m_eSd*sv.eSd(i,j)*sv.eSd(k,l); 
      }

      
      U3(i,j,k,l) = chi(i,j,k,l) + (-one_third*IoxI(i,j,r,s) - eSdoxeSd(i,j,r,s))*chi(r,s,k,l);
      
      // compute D(gamma_dot_v) part
      double factor3 = mat.param->gamma_dot_0/(mat.param->m*sv.g_pi)*pow_pi_pi_m_over_g_pi;

      Tensor<2> U4, U5, U6;
      U4(i,j) = dHdp*sv.pJ*one_third*I(k,l)*eFeSeFT(k,l)*sv.MI(j,i);
      U5(i,j) = -dHdp*sv.pJ*one_third*(2.0*Fa(k,i)*sv.eFnp1(k,l)*sv.eS(l,j) + eC(k,l)*chi(k,l,i,j));

      double dadpc = mat.compute_dadp(pc);
      double dcdp = mat.compute_dcdp(pc);      
      double dpi_mdp = dadpc - dcdp;
      double dg_pidp = mat.compute_b(sv.pi, sv.pi_m)*dadpc;
      
      U6(i,j) = (dpi_mdp + pi_pi_m_over_g_pi*dg_pidp)*sv.MI(j,i);
      
      // compute D(psi_v) part
      double factor4 = mat.param->g0*one_third/mat.param->pc_inf;
            
      U(i,j,k,l) = m_eSd*dHdp*A(i,k)*I(j,l) + dt*sv.eSd(i,j)*(U1(k,l) + U2(k,l))
                                 + dt*gamma_dot_d*U3(i,j,k,l)
                                 + dt*m_eSd*factor3*sv.psi_v(i,j)*(U4(k,l) + U5(k,l) + U6(k,l))
                                 - dt*m_eSd*gamma_dot_v*factor4*I(i,j)*sv.MI(l,k);
                                 
                                      
      // start computing B(right side);

      Tensor<4> U7;
      U7(i,j,k,l) = sv.L(i,j,k,l) + (-one_third*IoxI(i,j,r,s)-eSdoxeSd(i,j,r,s))*sv.L(r,s,k,l);
      
      for(int ia=0; ia<nne; ia++){
        for(int ib=0; ib<ndofn; ib++){
          int id_ab = ia*ndofn*DIM_3x3 + DIM_3x3*ib;
          
          TensorA<2> dMdu(dMdUs + id_ab), Grad_u(Grad_us + id_ab);
          
          Tensor<2> GradeFnMeSeF, GradeFnMeSeF_sym, eFGradeFnM, eFGradeFnM_sym;
          GradeFnMeSeF(i,j) = Grad_u(i,k)*eFn(k,l)*sv.M(l,m)*sv.eFnp1(j,m);
          eFGradeFnM(i,j) = sv.eFnp1(k,i)*Grad_u(k,l)*eFn(l,m)*sv.M(m,j);
          
          GradeFnMeSeF_sym(i,j) = 0.5*(GradeFnMeSeF(i,j) + GradeFnMeSeF(j,i));
          eFGradeFnM_sym(i,j) = 0.5*(eFGradeFnM(i,j) + eFGradeFnM(j,i));
          
          double factor5 = 2.0*d_bar_tau_d_tau(i,j)*GradeFnMeSeF_sym(i,j) + eFTdtaueF(i,j)*sv.L(i,j,k,l)*eFGradeFnM_sym(k,l);
          double factor6 = 2.0*I(i,j)*GradeFnMeSeF_sym(i,j) + eC(i,j)*sv.L(i,j,k,l)*eFGradeFnM_sym(k,l);
          B(i,j) = dHdp*(-dt*factor2*factor5*sv.eSd(i,j) - dt*gamma_dot_d*U7(i,j,k,l)*eFGradeFnM_sym(k,l)
                   + dt*m_eSd*sv.pJ*one_third*factor3*factor6*sv.psi_v(i,j));
          
          try{
            dMdu = ttl::solve( U, B);
          }
          catch(int solve_err){
            if(solver_info->debug)            
              printf("error on computing dMdu for pvp model. dMdu = delta_ij is set. Solution may not be converging.\n");
            dMdu(i,j) = O(i,j);
          }     
          
        }
      }

    } 

    ~PvpIntegrator(){};    
};

/// integration algorithm for PVP model in a staggered manner
///
/// \param[in]  param       material parameters(MaterialPoroViscoPlasticity) for PVP model 
/// \param[in]  solver_info contains numerical parameters such as number of maximum NR iteration and NR tolerance
/// \param[in]  Fnp1_in      F at t(n+1)
/// \param[in]  Fn_in        F at t(n)
/// \param[out] pFnp1_in    pF at t(n+1)
/// \param[in]  pFn_n       pF at t(n)
/// \param[out] pc_np1      conforming pressure at t(n+1)
/// \param[in]  pc_n        conforming pressure at t(n)
/// \param[in]  dt_in       time step size
/// \return non-zero with internal error
int poro_visco_plasticity_integration_algorithm_staggered(const MaterialPoroViscoPlasticity *param,
                                                          const GcmSolverInfo *solver_info,
                                                          double *Fnp1_in,
                                                          double *Fn_in,
                                                          double *pFnp1_in,
                                                          double *pFn_in,
                                                          double *pc_np1,
                                                          const double pc_n,
                                                          const double dt_in){
  int err = 0;
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(param); 
  
  try{      
    pvp.set_tenosrs(Fnp1_in, Fn_in, pFnp1_in, pFn_in); 
  }catch(int i){
    return 1;
  }
  
  pvp.set_scalars(*pc_np1, pc_n);
  pvp.set_solver_info(solver_info, dt_in);
                                                            
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
        if(solver_info->debug)
          printf("Matrix is singular. The solution (Poro-visco-plasticity) could not be computed.\n");
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

/// implicit integration algorithm for PVP model
///
/// \param[in]  param       material parameters(MaterialPoroViscoPlasticity) for PVP model 
/// \param[in]  solver_info contains numerical parameters such as number of maximum NR iteration and NR tolerance
/// \param[in]  Fnp1_in      F at t(n+1)
/// \param[in]  Fn_in        F at t(n)
/// \param[out] pFnp1_in    pF at t(n+1)
/// \param[in]  pFn_n       pF at t(n)
/// \param[out] pc_np1      conforming pressure at t(n+1)
/// \param[in]  pc_n        conforming pressure at t(n)
/// \param[in]  dt_in       time step size
/// \return non-zero with internal error
int poro_visco_plasticity_integration_algorithm_implicit(const MaterialPoroViscoPlasticity *param,
                                                         const GcmSolverInfo *solver_info,
                                                         double *Fnp1,
                                                         double *Fn,
                                                         double *pFnp1,
                                                         double *pFn,
                                                         double *pc_np1,
                                                         const double pc_n,
                                                         const double dt_in)
{
  int err = 0;
  
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(param);
  
  double dHdp = pvp.mat.compute_dHdp(pc_n);
  if(fabs(dHdp)<1.0e-6)
  {
    //printf("dHdp = %.17e\n", dHdp);  
    return poro_visco_plasticity_integration_algorithm_staggered(param, solver_info, Fnp1, Fn, pFnp1, pFn, pc_np1, pc_n, dt_in);
  }
  
  try{      
    pvp.set_tenosrs(Fnp1, Fn, pFnp1, pFn);
  }catch(int i){
    return 1;
  }
    
  pvp.set_scalars(*pc_np1, pc_n);
  pvp.set_solver_info(solver_info, dt_in);
  
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
      if(solver_info->debug)
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
      if(solver_info->debug)
        printf("converge on energe norm %e\n", eng_norm/eng_norm_0);

      //break;
    }
    
  } 
  *pc_np1 = pc;
  
  return err; 
}

/// exlicit integration algorithm for PVP model
///
/// \param[in]  param       material parameters(MaterialPoroViscoPlasticity) for PVP model 
/// \param[in]  solver_info contains numerical parameters such as number of maximum NR iteration and NR tolerance
/// \param[in]  Fnp1_in      F at t(n+1)
/// \param[in]  Fn_in        F at t(n)
/// \param[out] pFnp1_in    pF at t(n+1)
/// \param[in]  pFn_n       pF at t(n)
/// \param[out] pc_np1      conforming pressure at t(n+1)
/// \param[in]  pc_n        conforming pressure at t(n)
/// \param[in]  dt_in       time step size
/// \return non-zero with internal error
int poro_visco_plasticity_integration_algorithm_explicit(const MaterialPoroViscoPlasticity *param,
                                                         const GcmSolverInfo *solver_info,
                                                         double *Fnp1,
                                                         double *Fn,
                                                         double *pFnp1,
                                                         double *pFn,
                                                         double *pc_np1,
                                                         const double pc_n,
                                                         const double dt_in)
                                                         
{  
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(param);
  try{      
    pvp.set_tenosrs(Fnp1, Fn, pFnp1, pFn);
  }catch(int i){
    return 1;
  }
  pvp.set_scalars(*pc_np1, pc_n);
  pvp.set_solver_info(solver_info, dt_in);
  
  pvp.update_StateVariables(pc_n);
  Tensor<2> pD;
  TensorA<2> pF(pFnp1), pF_n(pFn);
  
  pvp.compute_pD(pD, pc_n);
  pF(i,j) = pF_n(i,j) + pvp.dt*pD(i,k)*pF_n(k,j);
  double pJ = ttl::det(pF);
  *pc_np1 = pvp.compute_pc(pc_n, pJ);
  
  return 0;
}

/// implicit integration algorithm for PVP model
///
/// \param[in]  param       material parameters(MaterialPoroViscoPlasticity) for PVP model 
/// \param[in]  solver_info contains numerical parameters such as number of maximum NR iteration and NR tolerance
/// \param[in]  Fnp1_in      F at t(n+1)
/// \param[in]  Fn_in        F at t(n)
/// \param[out] pFnp1_in    pF at t(n+1)
/// \param[in]  pFn_n       pF at t(n)
/// \param[out] pc_np1      conforming pressure at t(n+1)
/// \param[in]  pc_n        conforming pressure at t(n)
/// \param[in]  dt_in       time step size
/// \param[in]  is_implicit if yes: run implicit integration algorithm, default
/// \param[in]                 no : run explicit integration algorithm
/// \return non-zero with internal error
int poro_visco_plasticity_integration_algorithm(const MaterialPoroViscoPlasticity *param,
                                                const GcmSolverInfo *solver_info,
                                                double *Fnp1,
                                                double *Fn,
                                                double *pFnp1,
                                                double *pFn,
                                                double *pc_np1,
                                                const double pc_n,
                                                const double dt_in,
                                                const bool is_implicit)
{
  if(is_implicit)
    return poro_visco_plasticity_integration_algorithm_implicit(param, solver_info, Fnp1, Fn, pFnp1, pFn, pc_np1, pc_n, dt_in);
  else
    return poro_visco_plasticity_integration_algorithm_explicit(param, solver_info, Fnp1, Fn, pFnp1, pFn, pc_np1, pc_n, dt_in);
}

/// compute conforming pressure for a given plastic deformation 
///
/// \param[in] pJ        determinant of pF
/// \param[in] pc        initial comforming pressure
/// \param[in] param     poro_viscoplaticity material object
/// \return    computed  pc value  
double poro_visco_plasticity_compute_pc(double pJ, 
                                        double pc,
                                        const MaterialPoroViscoPlasticity *param,
                                        const GcmSolverInfo *solver_info)
{
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(param);
  pvp.set_solver_info(solver_info, 0.0);
        
  return pvp.compute_pc(pc, pJ);
}

/// compute PKII and elasticity tensor w.r.t eC
/// 
/// \param[out] eS_out computed PKII (2nd order tensor)
/// \param[out] L_out  computed elasticity tensor (4th order tensor)
/// \param[in]  param  poro_viscoplaticity material object  
/// \param[in]  eF_in  elastic part of deformation gradient
/// \param[in]  pc     conforming pressure
/// \param[in]  compute_4th_order if yes: compute L_out, default is false
///                                  no : skip computing L  
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

/// compute PKII and elasticity tensor for deviatoric part of W w.r.t eC
/// 
/// \param[out] eS_out computed PKII of hat_W (2nd order tensor)
/// \param[out] L_out  computed elasticity of tensor hat_W (4th order tensor)
/// \param[in]  param  poro_viscoplaticity material object  
/// \param[in]  eF_in  elastic part of deformation gradient
/// \param[in]  pc     conforming pressure
/// \param[in]  compute_4th_order if yes: compute L_out, default is false
///                                  no : skip computing L 
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

/// compute shear modulus as a function fo pc
/// 
/// \param[in]  param  poro_viscoplaticity material object  
/// \param[in]  pc     conforming pressure
/// \return computed shear modulus
double poro_visco_plasticity_compute_shear_modulus(const MaterialPoroViscoPlasticity *param,
                                                   const double pc){  
  PvpMaterial mat(param);      
  double c  = mat.compute_c(pc);
  double d  = mat.compute_d(pc);
  double mu = mat.compute_shear_modulus(c, d);
  
  return mu;
}

/// compute derivative of volumetric part of W(strain energy density function, U) w.r.t eJ
/// 
/// \param[in] eJ    det(eF)
/// \param[in] pc    conforming pressure
/// \param[in] param poro_viscoplaticity material object  
/// \return    computed dUdJ
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

/// compute 2nd derivative of volumetric part of W(strain energy density function, U) w.r.t eJ
/// 
/// \param[in] eJ    det(eF)
/// \param[in] pc    conforming pressure
/// \param[in] param poro_viscoplaticity material object  
/// \return    computed d2UdJ2
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

/// compute hardening function value of p_c 
///
/// \param[in] pc        initial comforming pressure
/// \param[in] param     poro_viscoplaticity material object
/// \return    computed H(p_c)
double poro_visco_plasticity_hardening(const double pc,
                                       const MaterialPoroViscoPlasticity *param){
  PvpMaterial mat(param);
  return mat.compute_H(pc);
}

/// compute plastic strain rates
///
/// \param[out] gamma_dot_d computed strain rate for deviatoric part
/// \param[out] gamma_dot_d computed strain rate for volumetric part
/// \param[in]   Fnp1_in     F at t(n+1)
/// \param[in]   Fn_in       F at t(n)
/// \param[in]  pFnp1_in    pF at t(n+1)
/// \param[in]  pFn_n       pF at t(n)
/// \param[in]  pc_np1      conforming pressure at t(n+1)
/// \param[in]  pc_n        conforming pressure at t(n)
/// \param[in]  dt_in       time step size
/// \param[in]  param       poro_viscoplaticity material object
/// \param[in]  solver_info contains numerical parameters such as number of maximum NR iteration and NR tolerance
void poro_visco_plasticity_intf_compute_gammas(double &gamma_dot_d,
                                               double &gamma_dot_v,
                                               double *Fnp1_in,
                                               double *Fn_in,
                                               double *pFnp1_in,
                                               double *pFn_in,
                                               const double pc_np1,
                                               const double pc_n,
                                               const double dt_in,
                                               const MaterialPoroViscoPlasticity *param,
                                               const GcmSolverInfo *solver_info)
{
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(param);       
  pvp.set_tenosrs(Fnp1_in, Fn_in, pFnp1_in, pFn_in); 

  pvp.set_scalars(pc_np1, pc_n);
  pvp.set_solver_info(solver_info, dt_in);
  pvp.update_StateVariables(pc_np1, false);
    
  double d = pvp.mat.compute_d(pc_np1);
  gamma_dot_d = pvp.mat.compute_gamma_dot_d(d, pvp.sv.bar_tau, pvp.sv.g_tau);
  gamma_dot_v = pvp.mat.compute_gamma_dot_v(pvp.sv.pi, pvp.sv.pi_m, pvp.sv.g_pi);  
}

/// compute dMdU for an finite element for PDE
/// 
/// \param[out] dMdUs       computed dMdU, nne*ndofn number of 2nd order tensors
/// \param[in]  Grad_us     Grad(delta u), nne*ndofn number of 2nd order tensors
/// \param[in]  param       poro_viscoplaticity material object
/// \param[in]  solver_info contains numerical parameters such as number of maximum NR iteration and NR tolerance
/// \param[in]   Fnp1_in     F at t(n+1)
/// \param[in]   Fn_in       F at t(n)
/// \param[in]  pFnp1_in    pF at t(n+1)
/// \param[in]  pFn_n       pF at t(n)
/// \param[in]  pc_np1      conforming pressure at t(n+1)
/// \param[in]  pc_n        conforming pressure at t(n)
/// \param[in]  dt          time step size
/// \param[in]  nne         number of nodes in an element
/// \param[in]  ndofn       number of nodal degree of freedoms
void poro_visco_plasticity_compute_dMdu(double *dMdUs,
                                        double *Grad_us,
                                        const MaterialPoroViscoPlasticity *param,
                                        const GcmSolverInfo *solver_info,
                                        double *Fnp1_in,
                                        double *Fn_in,
                                        double *pFnp1_in,
                                        double *pFn_in,
                                        const double pc_np1,
                                        const double pc_n,
                                        const double dt,
                                        const int nne,
                                        const int ndofn){
  PvpIntegrator pvp;
  pvp.set_pvp_material_parameters(param);       
  pvp.set_tenosrs(Fnp1_in, Fn_in, pFnp1_in, pFn_in, I.data, I.data, true);
  
  pvp.set_scalars(pc_np1, pc_n);
  pvp.set_solver_info(solver_info, dt);               
  
  pvp.update_StateVariables(pc_n, false);
  pvp.compute_dMdu(dMdUs, pc_np1, Grad_us, nne, ndofn);                                           
}

