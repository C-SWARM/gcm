#include "constitutive_model.h"
#include "continuum_damage_model.h"
#include "hyperelasticity.h"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))
  
constexpr double DELTA_W_MAX = 0.05;

int continuum_damage_Weibull(double *G, 
                             const double Y, 
                             const MATERIAL_CONTINUUM_DAMAGE *dam)
{
  int err = 0;  
  if(Y <= dam->Y_in) 
    *G = 0.0;
  else  
    *G = (dam->w_max - dam->w_max
	  * exp(- pow((Y - dam->Y_in)/(dam->P1*dam->Y_in),dam->P2)));  
  
  return err;
}

int continuum_damage_Weibull_evolution(double *H, 
                                       const double Y, 
                                       MATERIAL_CONTINUUM_DAMAGE const *dam)
{
  int err = 0;  
  if(Y <= dam->Y_in) 
    *H = 0.0;
  else  
    *H = (dam->w_max*dam->P2/(dam->P1*dam->Y_in)
	  *exp(-pow((Y - dam->Y_in) / (dam->P1 * dam->Y_in), dam->P2) )
	  * pow((Y - dam->Y_in)/(dam->P1 * dam->Y_in), dam->P2 - 1.0)
	  );  
  
  return err;
}

int damage_evolutions(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                      double *w,
                      double *X,
                      double *H,
                      int *is_it_damaged,
                      double wn,
                      double Xn,
                      const double Y,                                      
                      const double dt)
{
  int err = 0;
  const double dt_mu = dt*(mat_d->mu);
  double G = 0.0;
  
  err += continuum_damage_Weibull(&G,Y,mat_d);
  double g = G - Xn;

  if (g > 0.0)
  {    
    *is_it_damaged = 1;                      // Damage propagation
    *w = wn + dt_mu / (1 + dt_mu) * g;       // update damage parameter
    *X = MAX(Xn,(Xn + dt_mu*G)/(1.0+dt_mu)); // update softening parameter

    //update evolution parameter
    err += continuum_damage_Weibull_evolution(H,Y,mat_d);
  } 
  else 
  {
    *is_it_damaged = 0; //no damage propagation
    *w = wn;
    *X = Xn;
    *H = 0.0;
  }
  return err;
}

double H_of_J(double J)
{
  if(J>1.0)
    return 1.0;
  else
    return 0.0;
}
  
int split_damage_evolutions(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                            ELASTICITY *elast,
                            double *dw,
                            double *vw,
                            double *dX,
                            double *vX,
                            double *dH,
                            double *vH,
                            int *is_it_damaged_d,
                            int *is_it_damaged_v,                                     
                            double dwn,
                            double vwn,
                            double dXn,
                            double vXn,
                            const double W,
                            const double U,
                            const double J,
                            const double dt)
{
  int err = 0;
  const double dt_mu = dt*(mat_d->mu);
  double Gh = 0.0;
  double Gu = 0.0;
  
  double alpha_dev = mat_d->alpha_dev;
  double beta_dev  = mat_d->beta_dev;
  double alpha_vol = mat_d->alpha_vol;   
  double beta_vol  = mat_d->beta_vol;  
       
  if(beta_dev<0)
    beta_dev = -beta_dev*H_of_J(J);
  
  if(beta_vol<0)
    beta_vol = -beta_vol*H_of_J(J);      
          
  err += continuum_damage_Weibull(&Gh,alpha_dev*W + beta_dev*U,mat_d);
  err += continuum_damage_Weibull(&Gu,alpha_vol*W + beta_vol*U,mat_d);
    
  double gh = Gh - dXn;
  double gu = Gu - vXn;

  double H1, H2;
  err += continuum_damage_Weibull_evolution(&H1,alpha_dev*W + beta_dev*U,mat_d);
  err += continuum_damage_Weibull_evolution(&H2,alpha_vol*W + beta_vol*U,mat_d);

  
  if(gh > 0.0)
  {
    double dXnp1 = (dXn + dt_mu*Gh)/(1.0+dt_mu);
    *dH = H1*alpha_dev + H2*alpha_vol;
    
    *is_it_damaged_d = 1;
    *dw = dwn + dt_mu/(1.0+dt_mu)*(Gh - dXn);
    *dX = MAX(dXn,dXnp1);  
  }
  else 
  {
    *is_it_damaged_d = 0; //no damage propagation
    *dw = dwn;
    *dX = dXn;
    *dH = 0.0;
  }
  
  if(gu > 0.0)
  {
    double vXnp1 = (vXn + dt_mu*Gu)/(1.0+dt_mu);
    *vH = H1*beta_dev + H2*beta_vol;
    
    *is_it_damaged_v = 1;
    *vw = vwn + dt_mu/(1.0+dt_mu)*(Gu - vXn);
    *vX = MAX(vXn,vXnp1);    
  }
  else 
  {
    *is_it_damaged_v = 0; //no damage propagation
    *vw = vwn;
    *vX = vXn;
    *vH = 0.0;
  }
  return err;
}



int continuum_damage_integration_alg(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                     ELASTICITY *elast,
                                     double *w,
                                     double *X,
                                     double *H,
                                     int *is_it_damaged,
                                     double wn,
                                     double Xn,
                                     const double dt,
                                     double *F_in) 
{
  int err = 0;
  
  TensorA<2> F(F_in);
  double J = ttl::det(F);
    
  Tensor<2> C;
  C = F(k,i)*F(k,j);

  double Wdev = 0.0;
  double U = 0.0;
  
  
  
  elast->compute_potential_dev(C.data, elast->mat, &Wdev);
  elast->compute_u(&U, J);
  
  U *= (elast->mat->kappa);
  double Y = Wdev + U;
  

  err += damage_evolutions(mat_d,w,X,H,is_it_damaged,wn,Xn, Y, dt); 
  
  return err;
}

int continuum_damage_integration_alg_public(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                            ELASTICITY *elast,
                                            double *w,
                                            double *X,
                                            double *H,
                                            int *is_it_damaged,
                                            double wn,
                                            double Xn,
                                            const double dt,
                                            double Y)
{
  int err = 0;
      
  err += damage_evolutions(mat_d,w,X,H,is_it_damaged,wn,Xn, Y, dt);   
  return err;
}
int continuum_split_damage_integration_alg(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                           ELASTICITY *elast,
                                           double *dw,
                                           double *vw,
                                           double *dX,
                                           double *vX,
                                           double *dH,
                                           double *vH,
                                           int *is_it_damaged_d,
                                           int *is_it_damaged_v,                                     
                                           double dwn,
                                           double vwn,
                                           double dXn,
                                           double vXn,
                                           const double dt,
                                           double *F_in) 
{
  int err = 0;
  
  TensorA<2> F(F_in);
  double J = ttl::det(F);
    
  Tensor<2> C;
  C = F(k,i)*F(k,j);

  double W, U;
  W = U = 0.0;
  
  elast->compute_potential_dev(C.data, elast->mat, &W);    
  elast->compute_u(&U,J);
  U *= (elast->mat->kappa);
  
  err += split_damage_evolutions(mat_d,elast,dw,vw,dX,vX,dH,vH,
                                 is_it_damaged_d,is_it_damaged_v,
                                 dwn,vwn,dXn,vXn,
                                 W,U,J,dt);
  return err;
}

int continuum_split_damage_integration_alg_public(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                                  ELASTICITY *elast,
                                                  double *dw,
                                                  double *vw,
                                                  double *dX,
                                                  double *vX,
                                                  double *dH,
                                                  double *vH,
                                                  int *is_it_damaged_d,
                                                  int *is_it_damaged_v,                                     
                                                  double dwn,
                                                  double vwn,
                                                  double dXn,
                                                  double vXn,
                                                  const double dt,
                                                  double *F_in,
                                                  double W,
                                                  double U) 
{
  int err = 0;
  
  TensorA<2> F(F_in);
  double J = ttl::det(F);
    
  err += split_damage_evolutions(mat_d,elast,dw,vw,dX,vX,dH,vH,
                                 is_it_damaged_d,is_it_damaged_v,
                                 dwn,vwn,dXn,vXn,
                                 W,U,J,dt);    
  return err;
}

int apply_damage_on_stress(double *S, double *S0, double w)
{
  int err = 0;
  for(int a=0; a<DIM_3x3; a++)
    S[a] = (1.0 - w)*S0[a];
  return err;  
}

int apply_damage_on_stiffness(double *L_out, double *S0_in, double *L_in, 
                              double w, int is_it_damaged, double H,                               
                              double dt, double mu)
{
  int err = 0;
  for(int ia=0; ia<DIM_3x3x3x3; ia++)
    L_out[ia] *= (1.0 - w);
          
  if(is_it_damaged)
  {
    TensorA<2> S0(S0_in);
    TensorA<4> L(L_out);
    
    double dt_mu = dt*mu;
    double evo = dt_mu*H/(1.0+dt_mu);
    
    L(i,j,k,l) = L(i,j,k,l) - evo*S0(i,j)*S0(k,l);
  }            

  return err;  
}

int update_damage_elasticity(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                             ELASTICITY *elast,
                             double w,
                             int is_it_damaged,
                             double H,
                             const double dt,
                             double *F_in, 
                             const int compute_stiffness)
{
  int err = 0;
  elast->update_elasticity(elast,F_in, compute_stiffness);
              
  if(compute_stiffness)
  {
    TensorA<2> S(elast->S);
    Tensor<2> S0 = S;   
    
    for(int ia=0; ia<DIM_3x3x3x3; ia++)
      elast->L[ia] *= (1.0 - w);
         
    TensorA<4> L(elast->L);
    
    if(is_it_damaged)
    {
      double dt_mu = dt*(mat_d->mu);
      double evo = dt_mu*H/(1.0+dt_mu);
      
      L(i,j,k,l) = L(i,j,k,l) - evo*S0(i,j)*S0(k,l);
    }            
  }
  
  for(int a=0; a<DIM_3x3; a++)
    elast->S[a] *= (1.0 - w);
  
  return err;
}

int apply_split_damage_on_stress(double *S, double *dS0, double *vS0, double dw, double vw)
{
  int err = 0;
  for(int a=0; a<DIM_3x3; a++)
    S[a] = (1.0 - dw)*dS0[a] + (1.0 - vw)*vS0[a];
  return err;  
}

int apply_split_damage_on_stiffness(double *L_out, double *dS0_in, double *vS0_in,
                                    double *dL_in, double *vL_in, 
                                    double dw, double vw, int is_it_damaged_d, int is_it_damaged_v,
                                    double dH, double vH, double dt, double mu)
{
  int err = 0;
  
  TensorA<2> dS0(dS0_in), vS0(vS0_in);
  Tensor<4> dSS, vSS;
  
  dSS(i,j,k,l) = dS0(i,j)*dS0(k,l);
  vSS(i,j,k,l) = vS0(i,j)*vS0(k,l);
  
  double dt_mu = dt*mu;
      
  for(int a=0; a<DIM_3x3x3x3; a++)
  {    
    L_out[a] = (1-dw)*dL_in[a] + (1-vw)*vL_in[a];
    if(is_it_damaged_d)
      L_out[a] += (dt_mu)/(1.0+dt_mu)*(dH)*dSS.data[a];

    if(is_it_damaged_v)
      L_out[a] += (dt_mu)/(1.0+dt_mu)*(vH)*vSS.data[a];                           
  }
  return err;  
}

int update_split_damage_elasticity(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                   ELASTICITY *elast,
                                   double dw,
                                   double vw,
                                   double dH,
                                   double vH,
                                   int is_it_damaged_d,
                                   int is_it_damaged_v,
                                   const double dt,
                                   double *F_in, 
                                   const int compute_stiffness)
{
  int err = 0;

  Tensor<2> C, CI, dS_0;
  Tensor<2>vS_0 = {};
  
  // use double arrays as Matrix
  TensorA<2> F(F_in);
  C = F(k,i)*F(k,j);
  
  // compute stress -->
  elast->compute_PK2_dev(C.data, elast->mat, dS_0.data);

  inv(C,CI);
  double detF = ttl::det(F);
  double kappa = elast->mat->kappa;    
  double dudj = 0.0; 
  elast->compute_dudj(&dudj, detF);      

  vS_0(i,j) = kappa*detF*dudj*CI(i,j);  

  if(compute_stiffness)
  {            
    Tensor<4> dL;
    TensorA<4> L(elast->L);
    
    double detC = detF*detF;
    double d2udj2 = 0.0; 
    elast->compute_tangent_dev(C.data, elast->mat, dL.data);
    elast->compute_d2udj2(&d2udj2, detF);
    
    L(i,j,k,l) = (1-dw)*dL(i,j,k,l) + (1-vw)*(kappa*(detF*dudj + detC*d2udj2)*CI(i,j)*CI(k,l)
                                              -2.0*kappa*detF*dudj*CI(i,k)*CI(l,j));
 
    if(is_it_damaged_d || is_it_damaged_v)
    {  
      double W = 0.0, U = 0.0;

      elast->compute_potential_dev(C.data, elast->mat, &W);    
      elast->compute_u(&U,detF);
      U *= kappa;
    
      double dt_mu = dt*(mat_d->mu);
    
      if(is_it_damaged_d){
        double alpha_dev = mat_d->alpha_dev;
        double beta_dev  = mat_d->beta_dev;
      
        if(beta_dev<0)
          beta_dev = -beta_dev*H_of_J(detF);      
          
        double dY = alpha_dev*W + beta_dev*U;
        double H1 = 0.0;
        err += continuum_damage_Weibull_evolution(&H1, dY, mat_d);        
        L(i,j,k,l) += (dt_mu)/(1.0+dt_mu)*(H1*alpha_dev)*dS_0(i,j)*dS_0(k,l)
                   +  (dt_mu)/(1.0+dt_mu)*(H1*beta_dev)*dS_0(i,j)*vS_0(k,l);        
      }

      if(is_it_damaged_v){
        double alpha_vol = mat_d->alpha_vol;   
        double beta_vol  = mat_d->beta_vol;  
    
        if(beta_vol<0)
          beta_vol = -beta_vol*H_of_J(detF);

        double vY = alpha_vol*W + beta_vol*U;
        double H2 = 0.0;
        err += continuum_damage_Weibull_evolution(&H2, vY, mat_d);    
        L(i,j,k,l) += (dt_mu)/(1.0+dt_mu)*(H2*alpha_vol)*vS_0(i,j)*dS_0(k,l)
                   +  (dt_mu)/(1.0+dt_mu)*(H2*beta_vol)*vS_0(i,j)*vS_0(k,l);
      }
    }            
  }
  
  err += apply_split_damage_on_stress(elast->S, dS_0.data, vS_0.data, dw, vw);      
  return err;
}

int update_damage_elasticity_dev(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                 ELASTICITY *elast,
                                 double dw,
                                 double dH,
                                 int is_it_damaged_d,
                                 const double dt,
                                 double *F_in, 
                                 const int compute_stiffness)
{
  int err = 0;

  Tensor<2> C, dS_0;
  
  // use double arrays as Matrix
  TensorA<2> F(F_in);
  C = F(k,i)*F(k,j);
    
  // compute stress -->
  elast->compute_PK2_dev(C.data, elast->mat, dS_0.data);
          
  if(compute_stiffness)
  {
    Tensor<4> dL;
    TensorA<4> L(elast->L);
      
    elast->compute_tangent_dev(C.data, elast->mat, dL.data);              
    
    L(i,j,k,l) = (1-dw)*dL(i,j,k,l);
    
    if(is_it_damaged_d){
      Tensor<2> CI, vS_0;

      inv(C,CI);
      double detF = ttl::det(F);
      double kappa = elast->mat->kappa;    
      double dudj = 0.0; 
      elast->compute_dudj(&dudj, detF); 
  
      vS_0(i,j) = kappa*detF*dudj*CI(i,j);      
      
      double W = 0.0, U = 0.0;    
      elast->compute_potential_dev(C.data, elast->mat, &W);    
      elast->compute_u(&U,detF);
      U *= (elast->mat->kappa);
        
      double alpha_dev = mat_d->alpha_dev;
      double beta_dev  = mat_d->beta_dev;
           
      if(beta_dev<0)
        beta_dev = -beta_dev*H_of_J(detF);
      
        
      double dY = alpha_dev*W + beta_dev*U;  
      double H1 = 0.0;
      err += continuum_damage_Weibull_evolution(&H1, dY, mat_d);
      
      double dt_mu = dt*(mat_d->mu);
      L(i,j,k,l) += (dt_mu)/(1.0+dt_mu)*(H1*alpha_dev)*dS_0(i,j)*dS_0(k,l)
                 +  (dt_mu)/(1.0+dt_mu)*(H1*beta_dev)*dS_0(i,j)*vS_0(k,l);        
    }
  }  
  err += apply_damage_on_stress(elast->S, dS_0.data, dw); // apply only deviatoric part
  return err;
}


/// compute derivative of volumetric part of W(strain energy density function, U) w.r.t eJ
/// 
/// \param[in] elast elasticity object
/// \param[in] eJ    det(eF)
/// \param[in] vw    volumetric part damage parameter
double damage_compute_dudj(ELASTICITY *elast,
                           const double eJ,
                           const double vw)
{
  double dudj = 0.0;  
  elast->compute_dudj(&dudj, eJ);
  
  return (1.0 - vw)*dudj;  
}

/// compute  2nd derivative of volumetric part of W(strain energy density function, U) w.r.t eJ
/// 
/// \param[in] elast elasticity object
/// \param[in] eJ    det(eF)
/// \param[in] vw    volumetric part damage parameter
double damage_compute_d2udj2(ELASTICITY *elast,
                             const double eJ,
                             const double vw)
{
  double d2udj2 = 0.0;  
  elast->compute_d2udj2(&d2udj2, eJ);
  
  return (1.0 - vw)*d2udj2;  
}

/// Provide a public interface to the subdivision parameter to allow
/// other constitutive models to subdivide in a consistent fashion when
/// using the damage integration algorithm.
///
/// \param[in] w_n   damage at t(n)
/// \param[in] w_np1 damage at t(n+1)
double compute_subdiv_param(const double w_n,
                            const double w_np1){
  double w = w_np1 - w_n;
  return w/DELTA_W_MAX;
}


/// Provide a public interface to the subdivision parameter to allow
/// other constitutive models to subdivide in a consistent fashion when
/// using the split damage integration algorithm.
///
/// \param[in] dw_n deviatoric damage at t(n)
/// \param[in] dw_n deviatoric damage at t(n+1)
/// \param[in] vw_n volumetric damage at t(n)
/// \param[in] vw_n volumetric damage at t(n+1) 
double compute_subdiv_param4split_damage(const double dw_n,
                                         const double dw_np1,
                                         const double vw_n,
                                         const double vw_np1){
  double dw = dw_np1 - dw_n;
  double vw = vw_np1 - vw_n;
  return MAX(dw, vw)/DELTA_W_MAX;
}