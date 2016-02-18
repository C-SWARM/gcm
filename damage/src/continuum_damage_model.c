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

  dam->w  = w0;
  dam->X  = 0.0;
  dam->H  = 0.0;
  dam->wn = 0.0;
  dam->Xn = 0.0;
  dam->is_it_damaged = 0;
  
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
  const double dmu = dt*(dam->mat_d->mu);
  double G = 0.0;
  
  err += continuum_damage_Weibull(&G,Y,dam->mat_d);
  double g = G - dam->Xn;

  if (g > 0.0)
  {    
    dam->is_it_damaged = 1;                             //Damage propagation
    dam->w = dam->wn + dmu / (1 + dmu) * g;             // update damage parameter
    dam->X = MAX(dam->Xn,(dam->Xn + dmu*G)/ (1.0+dmu)); // update softening parameter

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
  
  double Y = Wdev + U*(elast->mat->kappa);

  damage_evolutions(dam,Y,dt);
  
  Matrix_cleanup(C); 
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

      double dmu = dt*(dam->mat_d->mu);
      double evo = dmu*(dam->H)/(1.0+dmu);
      
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