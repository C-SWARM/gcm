/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  Aaron Howell, [1], <ahowell3@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN

#include "constitutive_model.h"
#include "hyperelasticity.h"
#include "strain_energy_density_function.h"
#include "material_properties.h"
 
/// \param[out] elasticity
/// \param[in]  Fe
/// \param[in]  update_stiffness
int update_PK2_elasticity_tensor(ELASTICITY *elasticity, double *Fe, int update_stiffness)
{
  int err = 0;

  enum {C,CI,F2end};
  Tensor<2> F2[F2end];   //declare an array of tensors and initialize them to 0
  for (int a = 0; a < F2end; a++) {
    F2[a] = {};
  } 

  TensorA<2> F(Fe);
  TensorA<2> S(elasticity->S);

  F2[C](i,j) = F(k,i).to(i,k) * F(k,j);   //F[C] = F inverse * F
  
  err += inv(F2[C], F2[CI]);
    
  double detF = ttl::det(F);
  double detC = detF*detF;
  
  // compute stress -->
  double dudj = 0.0;  
  double kappa = elasticity->mat->kappa;    
  elasticity->compute_PK2_dev(F2[C].data, elasticity->mat, S.data);
  elasticity->compute_dudj(&dudj, detF);
  S(i,j) += kappa*detF*dudj * F2[CI](i,j);
  // <-- compute stress    

  //compute stiffness -->
  if(update_stiffness)
  {

    enum {CIoxCI,CICI,SoxS,F4end};
    Tensor<4> F4[F4end];   //declare an array of tensors and initialize them to 0
    for (int a = 0; a < F4end; a++) {
      F4[a] = {};
    }

    TensorA<4> L(elasticity->L);  
    
    double d2udj2 = 0.0;    
    elasticity->compute_tangent_dev(F2[C].data, elasticity->mat, L.data);
    elasticity->compute_d2udj2(&d2udj2, detF);

    F4[CIoxCI](i,j,k,l) = F2[CI](i,j) * F2[CI](k,l);   //calculate Kronecker product
    F4[SoxS]  (i,j,k,l) = S(i,j) * S(k,l);
    F4[CICI]  (i,j,k,l) = F2[CI](i,k) * F2[CI](l,j); 
  
    L(i,j,k,l) += kappa*(detF*dudj + detC*d2udj2) * F4[CIoxCI](i,j,k,l)
                - 2.0*kappa*detF*dudj * F4[CICI](i,j,k,l);   
  }
  // <-- compute stiffness
  
  return err;
}

int compute_effective_tensor2(double *T_in, double *T_eff)
{
  int err = 0;  
  
  TensorA<2> T(T_in);
  
  double seff = (T[0][0] - T[1][1])*(T[0][0] - T[1][1]) + 
                (T[1][1] - T[2][2])*(T[1][1] - T[2][2]) + 
                (T[2][2] - T[0][0])*(T[2][2] - T[0][0]) +
                6.0*(T[0][1]*T[0][1]+
                     T[1][2]*T[1][2]+
                     T[2][0]*T[2][0]);

  *T_eff = sqrt(seff/2.0);
  return err;
}

int compute_effective_PKII_stress(ELASTICITY *elasticity, double *PK2_eff)
{
  return compute_effective_tensor2(elasticity->S, PK2_eff);
}

int compute_Cauchy_stress(ELASTICITY *elasticity, double *sigma_out, double *eF)
{
  int err = 0;

  TensorA<2> Fe(eF), PK2(elasticity->S), sigma(sigma_out);        
  double det_Fe = ttl::det(Fe);  
  sigma = 1.0/det_Fe*Fe(i,k)*PK2(k,l)*Fe(j, l); 
    
  return err;  
}

int compute_effective_Cauchy_stress(ELASTICITY *elasticity, double *sigma_eff, double *eF)
{
  double *sigma = (double *) malloc(sizeof(double)*DIM_3x3);
  
  int err = compute_Cauchy_stress(elasticity, sigma, eF);
  err += compute_effective_tensor2(sigma, sigma_eff);
  return err;
}

int compute_tangent_of_tangent(ELASTICITY *elasticity, double *eF, double *K_in)
{
  int err = 0;
  
  TensorA<2> F(eF);
  Tensor<2> C, CI;  
  
  double J = ttl::det(F);
  C = F(k,i)*F(k,j);

  err += inv(C, CI);
  
  double kappa = elasticity->mat->kappa;    
  elasticity->compute_d3W_dC3_dev(C.data, elasticity->mat, K_in);
  TensorA<6> K(K_in);
  
  double dUdJ = 0.0;
  double d2UdJ2 = 0.0;
  double d3UdJ3 = 0.0;
  
  elasticity->compute_dudj(&dUdJ, J);
  elasticity->compute_d2udj2(&d2UdJ2, J);
  elasticity->compute_d3udj3(&d3UdJ3, J);
  
  K(i,j,k,l,r,s) += 2.0*((0.5*kappa*J*(dUdJ+3.0*J*d2UdJ2+J*J*d3UdJ3)*CI(i,j)*CI(k,l)*CI(r,s))
                          -(kappa*J*(dUdJ+J*d2UdJ2)*(CI(i,r)*CI(s,j)*CI(k,l) 
                                                   + CI(i,j)*CI(k,r)*CI(s,l) 
                                                   + CI(i,k)*CI(l,j)*CI(r,s)))
                          +2.0*kappa*J*dUdJ*(CI(i,r)*CI(s,k)*CI(l,j) + CI(i,k)*CI(l,r)*CI(s,j)));
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
