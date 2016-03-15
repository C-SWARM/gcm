#include "constitutive_model.h"
#include "continuum_damage_model.h"
#include "hyperelasticity.h"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

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
  
int damage_split_evolutions(MATERIAL_CONTINUUM_DAMAGE *mat_d,
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
  double g  = Gh + Gu;

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
  

  err += damage_evolutions(mat_d,w,X,H,is_it_damaged,wn,Xn, Y, dt); 
  
  Matrix_cleanup(C); 
  return err;
}

int continuum_damage_split_integration_alg(MATERIAL_CONTINUUM_DAMAGE *mat_d,
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
  
  Matrix(double) F;
  F.m_row = F.m_col = DIM_3; F.m_pdata = F_in;

  Matrix(double) C;
  Matrix_construct_redim(double,C    ,DIM_3,DIM_3);

  double W, U, J;
  W = U = J = 0.0;
  
  Matrix_AxB(C,1.0,0.0,F,1,F,0);
  Matrix_det(F, J);
  
  elast->compute_potential_dev(C.m_pdata, elast->mat, &W);    
  elast->compute_u(&U,J);
  U *= (elast->mat->kappa);
  
  err += damage_split_evolutions(mat_d,elast,dw,vw,dX,vX,dH,vH,
                                      is_it_damaged_d,is_it_damaged_v,
                                      dwn,vwn,dXn,vXn,
                                      W,U,J,dt);
    
  Matrix_cleanup(C);
  return err;
}

int update_damaged_elasticity(MATERIAL_CONTINUUM_DAMAGE *mat_d,
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

  Matrix(double) S;
  S.m_row = S.m_col = DIM_3;
  S.m_pdata = elast->S;
        
  if(compute_stiffness)
  {    
    Matrix(double) L;
    L.m_row = DIM_3x3x3x3; L.m_col = 1;
    L.m_pdata = elast->L;
    
    for(int I=1; I<=DIM_3x3x3x3; I++)
      Vec_v(L, I) *= (1.0 - w);
          
    if(is_it_damaged)
    {

      double dt_mu = dt*(mat_d->mu);
      double evo = dt_mu*H/(1.0+dt_mu);
      
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
    S.m_pdata[a] *= (1.0 - w);
  
  return err;
}


int update_damaged_elasticity_split(MATERIAL_CONTINUUM_DAMAGE *mat_d,
                                    ELASTICITY *elasticity,
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

  enum {C,CI,dS_0,vS_0,F2end};
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
  elasticity->compute_PK2_dev(F2[C].m_pdata, elasticity->mat, F2[dS_0].m_pdata);
  elasticity->compute_dudj(&dudj, detF);
  Matrix_AeqB(F2[vS_0], kappa*detF*dudj, F2[CI]);
          
  if(compute_stiffness)
  {
    enum {dL,dSS,vSS,CIoxCI,CICI,SoxS,F4end};
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
    elasticity->compute_tangent_dev(F2[C].m_pdata, elasticity->mat, F4[dL].m_pdata);              
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
            Tns4_v(F4[dSS]   ,I,J,P,Q) = Mat_v(F2[dS_0],I,J)*Mat_v(F2[dS_0],P,Q);
            Tns4_v(F4[vSS]   ,I,J,P,Q) = Mat_v(F2[vS_0],I,J)*Mat_v(F2[vS_0],P,Q);            
          }
        }
      }
    }
  
    double dt_mu = dt*(mat_d->mu);
        
    for(int I=1; I<=DIM_3x3x3x3; I++)
    {
      double vL = kappa*(detF*dudj + detC*d2udj2)*Vec_v(F4[CIoxCI], I)
                       - 2.0*kappa*detF*dudj*Vec_v(F4[CICI], I);
      Vec_v(L, I) = (1-dw)*Vec_v(F4[dL], I) + (1-vw)*vL;
      if(is_it_damaged_d)
        Vec_v(L, I) += (dt_mu)/(1.0+dt_mu)*(dH)*Vec_v(F4[dSS], I);

      if(is_it_damaged_v)
        Vec_v(L, I) += (dt_mu)/(1.0+dt_mu)*(vH)*Vec_v(F4[vSS], I);                           
    }

    for(int a = 0; a < F4end; a++)
      Matrix_cleanup(F4[a]);    
    free(F4);
  }
  
  for(int a=0; a<DIM_3x3; a++)
    elasticity->S[a] = (1.0 - dw)*F2[dS_0].m_pdata[a] 
                     + (1.0 - vw)*F2[vS_0].m_pdata[a];
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);    
  free(F2);
  
  return err;
}