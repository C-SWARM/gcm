/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  Aaron Howell, [1], <ahowell3@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN

#include "constitutive_model.h"
#include "hyperelasticity.h"
#include "strain_energy_density_function.h"
#include "material_properties.h"

constexpr const int Err = 1;
constexpr const double one_third = 1.0/3.0;
 
/// comute Von Mises stress
/// 
/// \param[in] T_in  2nd order tensor
/// \return computed Von_Mises value
double compute_Von_Mises(double *T_in)
{
  TensorA<2> T(T_in);
  
  double seff = (T[0][0] - T[1][1])*(T[0][0] - T[1][1]) + 
                (T[1][1] - T[2][2])*(T[1][1] - T[2][2]) + 
                (T[2][2] - T[0][0])*(T[2][2] - T[0][0]) +
                6.0*(T[0][1]*T[0][1]+
                     T[1][2]*T[1][2]+
                     T[2][0]*T[2][0]);

  return sqrt(seff/2.0);
}

/// set strain energy density function
/// \param[in] mat elastic material object
/// \return non-zero on internal error
int StrainEnergyDensityFunction::set_material(const MATERIAL_ELASTICITY *mat){
  int err = 0;

  switch (mat->devPotFlag)
  {
    case 1:
      this->compute_potential_dev = SEDF_devPotential_Mooney_Rivlin;
      this->compute_PK2_dev       = SEDF_devStress_Mooney_Rivlin;
      this->compute_tangent_dev   = SEDF_matStiff_Mooney_Rivlin;
      this->compute_d3W_dC3_dev   = SEDF_d3W_dC3_Mooney_Rivlin;
      break;
    case 2:
      this->compute_potential_dev = SEDF_devPotential_Linear;
      this->compute_PK2_dev       = SEDF_devStress_Linear;
      this->compute_tangent_dev   = SEDF_matStiff_Linear;
      this->compute_d3W_dC3_dev   = SEDF_d3W_dC3_Linear;      
      break;
    default:
      printf("ERROR: Unrecognized deviatoric potential flag (%d)\n", mat->devPotFlag);
      err = -1;
  }  
  
  switch(mat->volPotFlag)
  {
    case 99:
      this->compute_u      = SEDF_U_Common;
      this->compute_dudj   = SEDF_dUdJ_Common;
      this->compute_d2udj2 = SEDF_d2UdJ2_Common_new;
      this->compute_d3udj3 = SEDF_d3UdJ3_Common_new;
      break;  
    case 2:
      this->compute_u      = SEDF_U_Doll_Schweizerhof_7;
      this->compute_dudj   = SEDF_dUdJ_Doll_Schweizerhof_7;
      this->compute_d2udj2 = SEDF_d2UdJ2_Doll_Schweizerhof_7;
      this->compute_d3udj3 = SEDF_d3UdJ3_Doll_Schweizerhof_7;
      break;
    case 3:
      this->compute_u      = SEDF_U_Doll_Schweizerhof_8;
      this->compute_dudj   = SEDF_dUdJ_Doll_Schweizerhof_8;
      this->compute_d2udj2 = SEDF_d2UdJ2_Doll_Schweizerhof_8;
      this->compute_d3udj3 = SEDF_d3UdJ3_Doll_Schweizerhof_8;      
      break;
    default:
      printf("ERROR: Unrecognized volumetric potential flag (%d)\n",mat->volPotFlag);
      err = -1;
  }

  return err;
}

/// compute PKII stress and elasticity tensor (4th order) and update references (S_out, and L_out)
/// rather than member S and L 
/// if update_stiffness is false, there will be no computing elasticity tensor
/// s.t NULL is valid for L_out
/// 
/// \param[in]  mat   elastic material object
/// \param[out] S_out computed PKII
/// \param[out] L_out computed elasticity tensor
/// \param[in]  Fe    elastic deformation gradient
/// \param[in]  update_stiffness if true, compute 4th order elasticity Tensor
///                                 false no compute elasticity Tensor
/// \return non-zero on internal error
int StrainEnergyDensityFunction::update_elasticity(const MATERIAL_ELASTICITY *mat,
                                                   double *S_out,
                                                   double *L_out,
                                                   double *Fe, 
                                                   const bool update_stiffness){
  int err = 0;
  
  Tensor<2> C ={}, CI = {};

  TensorA<2> F(Fe);
  TensorA<2> eS(S_out);

  C(i,j) = F(k,i).to(i,k) * F(k,j);   //F[C] = F inverse * F
  
  err += inv(C, CI);
    
  double detF = ttl::det(F);
  double detC = detF*detF;
  
  // compute stress -->
  double dudj = 0.0;  
  double kappa = mat->kappa;    
  this->compute_PK2_dev(C.data, mat, eS.data);
  this->compute_dudj(&dudj, detF);
  eS(i,j) += kappa*detF*dudj * CI(i,j);
  // <-- compute stress    

  //compute stiffness -->
  if(update_stiffness)
  {
    Tensor<4> CIoxCI,CICI,SoxS,F4end;
    TensorA<4> LL(L_out);      
    double d2udj2 = 0.0;    
    this->compute_tangent_dev(C.data, mat, LL.data);
    this->compute_d2udj2(&d2udj2, detF);

    CIoxCI(i,j,k,l) = CI(i,j)*CI(k,l);   //calculate Kronecker product
    SoxS  (i,j,k,l) = eS(i,j)*eS(k,l);
    CICI  (i,j,k,l) = CI(i,k)*CI(l,j); 
  
    LL(i,j,k,l) += kappa*(detF*dudj + detC*d2udj2)*CIoxCI(i,j,k,l)
                - 2.0*kappa*detF*dudj*CICI(i,j,k,l);   
  }
  // <-- compute stiffness
  return err;
}

/// compute Cauchy stress from PKII in elastic configurateion
/// compute sigma = 1/eJ*eF*PKII*eF'
/// 
/// \param[in]  PKII      2nd Piolar Kirhhoff stress
/// \param[out] sigma_out computed Cauch stress
/// \param[in]  eF        elastic deformation tensor
void compute_Cauchy_stress(double *PKII,
                           double *sigma_out, 
                           double *eF){
  TensorA<2> Fe(eF), PK2(PKII), sigma(sigma_out);        
  double det_Fe = ttl::det(Fe);  
  sigma = 1.0/det_Fe*Fe(i,k)*PK2(k,l)*Fe(j, l);
}

/// compute 3rd derivative of strain energy function w.r.t eC
///
/// \param[in]  mat   elastic material object
/// \param[in]  eF    elastic deformation tensor
/// \param[out] K_out computed derivative 
/// \return non-zero on internal error of strain energy function w.r.t C
int StrainEnergyDensityFunction::compute_d3W_dC3(const MATERIAL_ELASTICITY *mat,
                                                 double *eF, 
                                                 double *K_out){  
  TensorA<2> F(eF);
  Tensor<2> C, CI;  
  
  double J = ttl::det(F);
  C = F(k,i)*F(k,j);

  int err = inv(C, CI);
  
  double kappa = mat->kappa;    
  this->compute_d3W_dC3_dev(C.data, mat, K_out);
  TensorA<6> K(K_out);
  
  double dUdJ = 0.0;
  double d2UdJ2 = 0.0;
  double d3UdJ3 = 0.0;
  
  this->compute_dudj(&dUdJ, J);
  this->compute_d2udj2(&d2UdJ2, J);
  this->compute_d3udj3(&d3UdJ3, J);
  
  K(i,j,k,l,r,s) += 2.0*((0.5*kappa*J*(dUdJ+3.0*J*d2UdJ2+J*J*d3UdJ3)*CI(i,j)*CI(k,l)*CI(r,s))
                          -(kappa*J*(dUdJ+J*d2UdJ2)*(CI(i,r)*CI(s,j)*CI(k,l) 
                                                   + CI(i,j)*CI(k,r)*CI(s,l) 
                                                   + CI(i,k)*CI(l,j)*CI(r,s)))
                          +2.0*kappa*J*dUdJ*(CI(i,r)*CI(s,k)*CI(l,j) + CI(i,k)*CI(l,r)*CI(s,j)));
  return err; 
}


