#include "constitutive_model.h"
#include "hyperelasticity.h"
#include "strain_energy_density_function.h"
#include "material_properties.h"
 
int update_PK2_elasticity_tensor(ELASTICITY *elasticity, double *Fe, int update_stiffness)
{
  int err = 0;

  // Matrix construct --->
  enum {C,CI,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_redim(double, F2[a],DIM_3,DIM_3);
  }
  // <-- Matrix construct
  
  // use double arrays as Matrix -->
  Matrix(double) F, S;
  F.m_row = F.m_col = DIM_3; F.m_pdata = Fe;
  S.m_row = S.m_col = DIM_3; S.m_pdata = elasticity->S;  
  // <-- use double array as Matrix
  
  double detF;
  Matrix_init(F2[C],0.0);
  Matrix_AxB(F2[C],1.0,0.0,F,1,F,0);  
  Matrix_inv(F2[C],F2[CI]);
  Matrix_det(F, detF);  
  double detC = detF*detF;
  
  // compute stress -->
  double dudj = 0.0;  
  double kappa = elasticity->mat->kappa;    
  elasticity->compute_PK2_dev(F2[C].m_pdata, elasticity->mat, S.m_pdata);
  elasticity->compute_dudj(&dudj, detF);
  Matrix_AplusB(S, kappa*detF*dudj,F2[CI],1.0,S);
  // <-- compute stress    

  //compute stiffness -->
  if(update_stiffness)
  {
    
    enum {CIoxCI,CICI,SoxS,F4end};
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
    elasticity->compute_tangent_dev(F2[C].m_pdata, elasticity->mat, L.m_pdata);              
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
          }
        }
      }
    }
  
    for(int I=1; I<=DIM_3x3x3x3; I++)
    {
      Vec_v(L, I) += kappa*(detF*dudj + detC*d2udj2)*Vec_v(F4[CIoxCI], I)
                   - 2.0*kappa*detF*dudj*Vec_v(F4[CICI], I);
    }    

    for(int a = 0; a < F4end; a++)
      Matrix_cleanup(F4[a]);    
    free(F4);
  }
  // <-- compute stiffness
    
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);    
  free(F2);
  
  return 0;
}

int compute_effective_tensor2(double *T_in, double *T_eff)
{
  int err = 0;  
  
  Matrix(double) T;
  T.m_row = T.m_col = DIM_3; T.m_pdata = T_in; 
  
  double seff = (Mat_v(T,1,1) - Mat_v(T,2,2))*(Mat_v(T,1,1) - Mat_v(T,2,2)) + 
                (Mat_v(T,2,2) - Mat_v(T,3,3))*(Mat_v(T,2,2) - Mat_v(T,3,3)) + 
                (Mat_v(T,3,3) - Mat_v(T,1,1))*(Mat_v(T,3,3) - Mat_v(T,1,1)) +
                6.0*(Mat_v(T,1,2)*Mat_v(T,1,2)+
                     Mat_v(T,2,3)*Mat_v(T,2,3)+
                     Mat_v(T,3,1)*Mat_v(T,3,1));

  *T_eff = sqrt(seff/2.0);
  return err;

/*                
  double trT;
  Matrix_trace(T,trT);

  Matrix(double) Tdev;
  Matrix_construct_redim(double,Tdev,DIM_3,DIM_3);
  
  Matrix_eye(Tdev, 3);                
  Matrix_AplusB(Tdev, 1.0, T, -trT/3.0, Tdev);
    
  double norm_T;
  Matrix_ddot(Tdev,Tdev,norm_T);
  *T_eff = sqrt(3.0/2.0*norm_T);

  Matrix_cleanup(Tdev); 

  return err;    */
}

int compute_effective_PKII_stress(ELASTICITY *elasticity, double *PK2_eff)
{
  return compute_effective_tensor2(elasticity->S, PK2_eff);
}

int compute_Cauchy_stress(ELASTICITY *elasticity, double *sigma_out, double *eF)
{
  int err = 0;

  Matrix(double) PK2, Fe, sigma;
     Fe.m_row =    Fe.m_col = DIM_3;    Fe.m_pdata = eF;   
    PK2.m_row =   PK2.m_col = DIM_3;   PK2.m_pdata = elasticity->S; 
  sigma.m_row = sigma.m_col = DIM_3; sigma.m_pdata = sigma_out; 
    
  Matrix(double) FeT;
  Matrix_construct_redim(double,  FeT,DIM_3,DIM_3);
  
  Matrix_AeqBT(FeT,1.0, Fe);   
 
  double det_Fe;
  Matrix_det(Fe, det_Fe);
   
  Matrix_Tns2_AxBxC(sigma,1.0/det_Fe,0.0,Fe,PK2,FeT);
  Matrix_cleanup(FeT);
    
  return err;  
}

int compute_effective_Cauchy_stress(ELASTICITY *elasticity, double *sigma_eff, double *eF)
{
  double *sigma = (double *) malloc(sizeof(double)*DIM_3x3);
  
  int err = compute_Cauchy_stress(elasticity, sigma, eF);
  err += compute_effective_tensor2(sigma, sigma_eff);
  return err;
}

int compute_tangent_of_tangent(ELASTICITY *elasticity, double *eF, double *K)
{
  int err = 0;

  // Matrix construct --->
  enum {C,CI,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_redim(double, F2[a],DIM_3,DIM_3);
  }
  // <-- Matrix construct
  
  Matrix(double) F;
  F.m_row = F.m_col = DIM_3; F.m_pdata = eF;
  
  double J = 0.0;
  Matrix_init(F2[C],0.0);
  Matrix_AxB(F2[C],1.0,0.0,F,1,F,0);  
  Matrix_inv(F2[C],F2[CI]);
  Matrix_det(F, J); 
  
  double *C_I = F2[CI].m_pdata; 
  
  double kappa = elasticity->mat->kappa;    
  elasticity->compute_d3W_dC3_dev(F2[C].m_pdata, elasticity->mat, K);
  
  double dUdJ = 0.0;
  double d2UdJ2 = 0.0;
  double d3UdJ3 = 0.0;
  
  elasticity->compute_dudj(&dUdJ, J);
  elasticity->compute_d2udj2(&d2UdJ2, J);
  elasticity->compute_d3udj3(&d3UdJ3, J);
  
  for(int i=0; i<DIM_3; i++)
  {
    for(int j=0; j<DIM_3; j++)
    {
      const int ij = DIM_3*i+j;
      for(int k=0; k<DIM_3; k++)
      {
        const int ik = DIM_3*i+k;
        for(int l=0; l<DIM_3; l++)
        {
          const int kl = DIM_3*k + l;
          const int lj = DIM_3*l + j;
          const int ijkl = DIM_3x3x3*i+DIM_3x3*j + DIM_3*k+l;
          
          for(int r=0; r<DIM_3; r++)
          {
            const int ir = DIM_3*i + r;
            const int kr = DIM_3*k + r;
            const int lr = DIM_3*l + r;
            
            for(int s=0; s<DIM_3; s++)
            {
              const int sj = DIM_3*s + j;
              const int sk = DIM_3*s + k;
              const int sl = DIM_3*s + l;
              const int rs = DIM_3*r + s;
              const int ijrs = DIM_3x3x3*i + DIM_3x3*j + DIM_3*r + s;
              const int klrs = DIM_3x3x3*k + DIM_3x3*l + DIM_3*r + s;
              const int ijklrs = DIM_3x3x3x3*3*i + DIM_3x3x3x3*j + DIM_3x3x3*k + DIM_3x3*l + DIM_3*r + s;
              
              const double K_vol = 2.*(
                      (0.5*kappa*J*(dUdJ+3*J*d2UdJ2+J*J*d3UdJ3)
                      *C_I[ij]*C_I[kl]*C_I[rs])
                      
                      -(kappa*J*(dUdJ+J*d2UdJ2)
                      *(C_I[ir]*C_I[sj]*C_I[kl]
                      + C_I[ij]*C_I[kr]*C_I[sl]
                      + C_I[ik]*C_I[lj]*C_I[rs]))
                      
                      +2.*kappa*J*dUdJ*(C_I[ir]*C_I[sk]*C_I[lj]
                      + C_I[ik]*C_I[lr]*C_I[sj])
                      );
              
              K[ijklrs] += K[ijklrs];
              
            }
          }
        }
      }
    }
  }
    
  return err;  
}

int construct_elasticity(ELASTICITY *elasticity, MATERIAL_ELASTICITY *mat, const int compute_stiffness)
{
  int err = 0;  
  elasticity->S = (double *) malloc(sizeof(double)*DIM_3x3);
  if(compute_stiffness)  
    elasticity->L = (double *) malloc(sizeof(double)*DIM_3x3x3x3);
  elasticity->update_elasticity = update_PK2_elasticity_tensor;
  elasticity->mat = mat;
  elasticity->compute_stiffness = compute_stiffness;

  elasticity->compute_PK2_eff     = compute_effective_PKII_stress;
  elasticity->compute_Cauchy_eff  = compute_effective_Cauchy_stress;
  elasticity->compute_Cauchy      = compute_Cauchy_stress; 
  elasticity->compute_d3W_dC3     = compute_tangent_of_tangent;  
  switch (mat->devPotFlag)
  {
    case 1:
      elasticity->compute_potential_dev = SEDF_devPotential_Mooney_Rivlin;
      elasticity->compute_PK2_dev     = SEDF_devStress_Mooney_Rivlin;
      elasticity->compute_tangent_dev = SEDF_matStiff_Mooney_Rivlin;
      elasticity->compute_d3W_dC3_dev = SEDF_d3W_dC3_Mooney_Rivlin;
      break;
    case 2:
      elasticity->compute_potential_dev = SEDF_devPotential_Linear;
      elasticity->compute_PK2_dev     = SEDF_devStress_Linear;
      elasticity->compute_tangent_dev = SEDF_matStiff_Linear;
      elasticity->compute_d3W_dC3_dev = SEDF_d3W_dC3_Linear;      
      break;
    default:
      printf("ERROR: Unrecognized deviatoric potential flag (%d)\n", mat->devPotFlag);
      err = -1;
  }  
  
  switch(mat->volPotFlag)
  {
    case 99:
      elasticity->compute_u      = SEDF_U_Common;
      elasticity->compute_dudj   = SEDF_dUdJ_Common;
      elasticity->compute_d2udj2 = SEDF_d2UdJ2_Common_new;
      elasticity->compute_d3udj3 = SEDF_d3UdJ3_Common_new;
      break;  
    case 2:
      elasticity->compute_u      = SEDF_U_Doll_Schweizerhof_7;
      elasticity->compute_dudj   = SEDF_dUdJ_Doll_Schweizerhof_7;
      elasticity->compute_d2udj2 = SEDF_d2UdJ2_Doll_Schweizerhof_7;
      elasticity->compute_d3udj3 = SEDF_d3UdJ3_Doll_Schweizerhof_7;
      break;
    case 3:
      elasticity->compute_u      = SEDF_U_Doll_Schweizerhof_8;
      elasticity->compute_dudj   = SEDF_dUdJ_Doll_Schweizerhof_8;
      elasticity->compute_d2udj2 = SEDF_d2UdJ2_Doll_Schweizerhof_8;
      elasticity->compute_d3udj3 = SEDF_d3UdJ3_Doll_Schweizerhof_8;      
      break;
    default:
      printf("ERROR: Unrecognized volumetric potential flag (%d)\n",mat->volPotFlag);
      err = -1;
  }
    
  return err;
}

int destruct_elasticity(ELASTICITY *elasticity)
{
  int err = 0; 
  elasticity->mat                 = NULL; 
  elasticity->update_elasticity   = NULL;
  elasticity->compute_d3W_dC3     = NULL;
  elasticity->compute_PK2_dev     = NULL;
  elasticity->compute_tangent_dev = NULL;
  elasticity->compute_dudj        = NULL;
  elasticity->compute_d2udj2      = NULL;
  elasticity->compute_PK2_eff     = NULL;
  elasticity->compute_Cauchy_eff  = NULL;
  elasticity->compute_Cauchy      = NULL;   
  
  free(elasticity->S);
  if(elasticity->compute_stiffness)
    free(elasticity->L);
    
  return err;
}
