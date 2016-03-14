#include "constitutive_model.h"
#include "continuum_damage_model.h"
#include "hyperelasticity.h"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

int initialize_continuum_damage(CONTINUUM_DAMAGE *dam, 
                                ELASTICITY *elast, 
                                MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                double w0)
{
  int err = 0;

  dam->w  = 0.0;
  dam->X  = 0.0;
  dam->H  = 0.0;
  dam->wn = 0.0;
  dam->Xn = w0;
  dam->is_it_damaged = 0;
  
  dam->elasticity = elast;
  dam->mat_d      = mat_d; 
  return err;
}
                                
int initialize_continuum_damage_split(CONTINUUM_DAMAGE_SPLIT *dam, 
                                      ELASTICITY *elast, 
                                      MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                      double w0)
{
  int err = 0;

  dam->dw   = 0.0;
  dam->vw   = 0.0;
  dam->dX   = 0.0;
  dam->vX   = 0.0;    
  dam->dH   = 0.0;
  dam->vH   = 0.0;
  dam->dwn = w0;  
  dam->vwn = w0;
  dam->dXn = 0.0;
  dam->vXn = 0.0;    
  dam->is_it_damaged_h = 0;
  dam->is_it_damaged_u = 0;
  
  dam->elasticity = elast;
  dam->mat_d      = mat_d; 
  return err;  
}                                      
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

int damage_evolutions(CONTINUUM_DAMAGE *dam,
                      const double Y,                                      
                      const double dt)
{
  int err = 0;
  const double dt_mu = dt*(dam->mat_d->mu);
  double G = 0.0;
  
  err += continuum_damage_Weibull(&G,Y,dam->mat_d);
  double g = G - dam->Xn;

  if (g > 0.0)
  {    
    dam->is_it_damaged = 1;                             //Damage propagation
    dam->w = dam->wn + dt_mu / (1 + dt_mu) * g;             // update damage parameter
    dam->X = MAX(dam->Xn,(dam->Xn + dt_mu*G)/(1.0+dt_mu)); // update softening parameter

    //update evolution parameter
    err += continuum_damage_Weibull_evolution(&(dam->H),Y,dam->mat_d);
  } 
  else 
  {
    dam->is_it_damaged = 0; //no damage propagation
    dam->w = dam->wn;
    dam->X = dam->Xn;
    dam->H = 0.0;
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
  
int damage_split_evolutions(CONTINUUM_DAMAGE_SPLIT *dam,
                            const double W,
                            const double U,
                            const double J,
                            const double dt)
{
  int err = 0;
  const double dt_mu = dt*(dam->mat_d->mu);
  double Gh = 0.0;
  double Gu = 0.0;
   
  double alpha_dev = (dam->mat_d)->alpha_dev;
  double beta_dev  = (dam->mat_d)->beta_dev;
  double alpha_vol = (dam->mat_d)->alpha_vol;   
  double beta_vol  = (dam->mat_d)->beta_vol;
    
  if(beta_dev<0)
    beta_dev = -beta_dev*H_of_J(J);
  
  if(beta_vol<0)
    beta_vol = -beta_vol*H_of_J(J);      
          
  err += continuum_damage_Weibull(&Gh,alpha_dev*W + beta_dev*U,dam->mat_d);
  err += continuum_damage_Weibull(&Gu,alpha_vol*W + beta_vol*U,dam->mat_d);
    
  double gh = Gh - dam->dXn;
  double gu = Gu - dam->vXn;
  double g  = Gh + Gu;

  double H1, H2;
  err += continuum_damage_Weibull_evolution(&H1,alpha_dev*W + beta_dev*U,dam->mat_d);
  err += continuum_damage_Weibull_evolution(&H2,alpha_vol*W + beta_vol*U,dam->mat_d);

  
  if(gh > 0.0)
  {
    double dX = (dam->dXn + dt_mu*Gh)/(1.0+dt_mu);
    double dH = H1*alpha_dev + H2*alpha_vol;
    
    dam->is_it_damaged_h = 1;
    dam->dw = dam->dwn + dt_mu/(1.0+dt_mu)*(Gh - dam->dXn);
    dam->dX = MAX(dam->dXn,dX);  
  }
  else 
  {
    dam->is_it_damaged_h = 0; //no damage propagation
    dam->dw = dam->dwn;
    dam->dX = dam->dXn;
    dam->dH = 0.0;
  }
  
  if(gu > 0.0)
  {
    double vX = (dam->vXn + dt_mu*Gu)/(1.0+dt_mu);
    double vH = H1*beta_dev + H2*beta_vol;
    
    dam->is_it_damaged_u = 1;
    dam->vw = dam->vwn + dt_mu/(1.0+dt_mu)*(Gu - dam->vXn);
    dam->vX = MAX(dam->vXn,vX);  
  }
  else 
  {
    dam->is_it_damaged_u = 0; //no damage propagation
    dam->vw = dam->vwn;
    dam->vX = dam->vXn;
    dam->vH = 0.0;
  }
/*    
//    if(J>1.0)
//    {

      err += continuum_damage_Weibull_evolution(&vH,U,dam->mat_d);      

      dam->is_it_damaged_u = 1;
//    }
//    else
//    {
//      err += continuum_damage_Weibull_evolution(&H,W,dam->mat_d);
//      dH = H;
//      vH = 0.0;
//    }
    
    dam->dw = dam->dwn + dt_mu*dH*(Gh+Gu-X)/(dH+vH);
    dam->vw = dam->vwn + dt_mu*vH*(Gh+Gu-X)/(dH+vH);    
    
    dam->X = MAX(dam->Xn,X); // update softening parameter    
    dam->dH = dH;
    dam->vH = vH;

  } 
  else 
  {
    dam->is_it_damaged_h = 0; //no damage propagation
    dam->is_it_damaged_u = 0; //no damage propagation    
    dam->dw = dam->dwn;
    dam->vw = dam->vwn;
    dam->X = dam->Xn;
    dam->dH = 0.0;
    dam->vH = 0.0;
  } */
  return err;
}

int continuum_damage_integration_alg(CONTINUUM_DAMAGE *dam,
                                     const double dt,
                                     double *F_in) 
{
  int err = 0;
  
  ELASTICITY *elast = dam->elasticity;
  
  Matrix(double) F;
  F.m_row = F.m_col = DIM_3; F.m_pdata = F_in;

  Matrix(double) C;
  Matrix_construct_redim(double,C,DIM_3,DIM_3);

  double Wdev = 0.0;
  double U = 0.0;
  
  Matrix_AxB(C,1.0,0.0,F,1,F,0);
  double J = 0.0;
  Matrix_det(F,J);
  
  elast->compute_potential_dev(C.m_pdata, elast->mat, &Wdev);
  elast->compute_u(&U, J);
  
  U *= (elast->mat->kappa);
  double Y = Wdev + U;
  
  damage_evolutions(dam,Y,dt);
  
  Matrix_cleanup(C); 
  return err;
}

int continuum_damage_split_integration_alg(CONTINUUM_DAMAGE_SPLIT *dam,
                                           const double dt,
                                           double *F_in) 
{
  int err = 0;
  
  ELASTICITY *elast = dam->elasticity;
  
  Matrix(double) F;
  F.m_row = F.m_col = DIM_3; F.m_pdata = F_in;

  Matrix(double) C, sigma;
  Matrix_construct_redim(double,C    ,DIM_3,DIM_3);
  Matrix_construct_redim(double,sigma,DIM_3,DIM_3);  

  double W, U, J;
  W = U = J = 0.0;
  
  Matrix_AxB(C,1.0,0.0,F,1,F,0);
  Matrix_det(F, J);
  
  elast->compute_potential_dev(C.m_pdata, elast->mat, &W);    
  elast->compute_u(&U,J);
  U *= (elast->mat->kappa);
  
  elast->compute_Cauchy(elast,sigma.m_pdata,F_in);
  damage_split_evolutions(dam, W, U, J, dt);
  
  Matrix_cleanup(C);
  Matrix_cleanup(sigma);
  return err;
}

int update_damage_time_steps(CONTINUUM_DAMAGE *dam)
{
  int err = 0;
  dam->wn = dam->w;
  dam->Xn = dam->X;
  dam->is_it_damaged = 0;    
  return err;
}

int update_damage_split_time_steps(CONTINUUM_DAMAGE_SPLIT *dam)
{
  int err = 0;
  dam->dwn = dam->dw;
  dam->vwn = dam->vw;   
  dam->dXn = dam->dX;
  dam->vXn = dam->vX;  
  dam->is_it_damaged_h = 0;    
  dam->is_it_damaged_u = 0;  
  return err;
}

int update_damaged_elasticity(CONTINUUM_DAMAGE *dam, 
                              const double dt,
                              double *F_in, 
                              const int compute_stiffness)
{
  int err = 0;
  ELASTICITY *elast = dam->elasticity;
  elast->update_elasticity(elast,F_in, compute_stiffness);

  Matrix(double) S;
  S.m_row = S.m_col = DIM_3;
  S.m_pdata = elast->S;
        
  if(compute_stiffness)
  {    
    Matrix(double) L;
    L.m_row = DIM_3x3x3x3; L.m_col = 1;
    L.m_pdata = elast->L;
    
    for(int I=1; I<=DIM_3x3x3x3; I++)
      Vec_v(L, I) *= (1.0 - dam->w);
          
    if(dam->is_it_damaged)
    {

      double dt_mu = dt*(dam->mat_d->mu);
      double evo = dt_mu*(dam->H)/(1.0+dt_mu);
      
      Matrix(double) S0;
      Matrix_construct_redim(double, S0,DIM_3,DIM_3);
      Matrix_AeqB(S0,1.0,S);
      
      for(int I=1; I<=DIM_3; I++)
        for(int J=1; J<=DIM_3; J++)
          for(int P=1; P<=DIM_3; P++)
            for(int Q=1; Q<=DIM_3; Q++)
              Tns4_v(L,I,J,P,Q) -= evo*Mat_v(S0,I,J)*Mat_v(S0,P,Q);     
      
      Matrix_cleanup(S0);
    }            
  }
  
  for(int a=0; a<DIM_3x3; a++)
    S.m_pdata[a] *= (1.0 - dam->w);
  
  return err;
}

int update_damaged_elasticity_split(CONTINUUM_DAMAGE_SPLIT *dam, 
                              const double dt,
                              double *F_in, 
                              const int compute_stiffness)
{
  int err = 0;
  ELASTICITY *elasticity = dam->elasticity;

  enum {C,CI,Sh_0,Su_0,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_redim(double, F2[a],DIM_3,DIM_3);
  }
  // <-- Matrix construct
  
  // use double arrays as Matrix -->
  Matrix(double) F, S;
  F.m_row = F.m_col = DIM_3; F.m_pdata = F_in;
  S.m_row = S.m_col = DIM_3; S.m_pdata = elasticity->S;  
  // <-- use double array as Matrix
  
  double detF;
  Matrix_AxB(F2[C],1.0,0.0,F,1,F,0);  
  Matrix_inv(F2[C],F2[CI]);
  Matrix_det(F, detF);  
  double detC = detF*detF;
  
  // compute stress -->
  double dudj = 0.0;  
  double kappa = elasticity->mat->kappa;    
  elasticity->compute_PK2_dev(F2[C].m_pdata, elasticity->mat, F2[Sh_0].m_pdata);
  elasticity->compute_dudj(&dudj, detF);
  Matrix_AeqB(F2[Su_0], kappa*detF*dudj, F2[CI]);
          
  if(compute_stiffness)
  {
    enum {Lh,SSh,SSu,CIoxCI,CICI,SoxS,F4end};
    Matrix(double) *F4 = malloc(F4end*sizeof(Matrix(double)));
    for (int a = 0; a < F4end; a++) {
      Matrix_construct_redim(double, F4[a],DIM_3x3x3x3,1);
    }
    // use double arrays as Matrix -->
    Matrix(double) L;
    L.m_row = DIM_3x3x3x3; L.m_col = 1;
    L.m_pdata = elasticity->L;
      
    // <-- use double array as Matrix        
    
    double d2udj2 = 0.0;    
    elasticity->compute_tangent_dev(F2[C].m_pdata, elasticity->mat, F4[Lh].m_pdata);              
    elasticity->compute_d2udj2(&d2udj2, detF);

    for(int I=1; I<=DIM_3; I++)
    {
      for(int J=1; J<=DIM_3; J++)
      {
        for(int P=1; P<=DIM_3; P++)
        {
          for(int Q=1; Q<=DIM_3; Q++)
          {
            Tns4_v(F4[CIoxCI],I,J,P,Q) = Mat_v(F2[CI],I,J)*Mat_v(F2[CI],P,Q);
            Tns4_v(F4[SoxS]  ,I,J,P,Q) = Mat_v(S,I,J)*Mat_v(S,P,Q);            
            Tns4_v(F4[CICI]  ,I,J,P,Q) = Mat_v(F2[CI],I,P)*Mat_v(F2[CI],Q,J);
            Tns4_v(F4[SSh]   ,I,J,P,Q) = Mat_v(F2[Sh_0],I,J)*Mat_v(F2[Sh_0],P,Q);
            Tns4_v(F4[SSu]   ,I,J,P,Q) = Mat_v(F2[Su_0],I,J)*Mat_v(F2[Su_0],P,Q);            
          }
        }
      }
    }
  
    double dt_mu = dt*(dam->mat_d->mu);
        
    for(int I=1; I<=DIM_3x3x3x3; I++)
    {
      double Lu = kappa*(detF*dudj + detC*d2udj2)*Vec_v(F4[CIoxCI], I)
                       - 2.0*kappa*detF*dudj*Vec_v(F4[CICI], I);
      Vec_v(L, I) = (1-dam->dw)*Vec_v(F4[Lh], I) + (1-dam->dw)*Lu;
      if(dam->is_it_damaged_h)
        Vec_v(L, I) += (dt_mu)/(1.0+dt_mu)*(dam->dH)*Vec_v(F4[SSh], I);

      if(dam->is_it_damaged_u)
        Vec_v(L, I) += (dt_mu)/(1.0+dt_mu)*(dam->vH)*Vec_v(F4[SSu], I);                           
    }

    for(int a = 0; a < F4end; a++)
      Matrix_cleanup(F4[a]);    
    free(F4);
  }
  
  for(int a=0; a<DIM_3x3; a++)
    elasticity->S[a] = (1.0 - dam->dw)*F2[Sh_0].m_pdata[a] 
                     + (1.0 - dam->vw)*F2[Su_0].m_pdata[a];
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);    
  free(F2);
  
  return err;
}