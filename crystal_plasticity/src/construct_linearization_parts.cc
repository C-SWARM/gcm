/// Authors:
///  Sangmin Lee, [1], <slee43@nd.edu>
///  Aaron Howell, [1], <ahowell3@nd.edu>
///  [1] - University of Notre Dame, Notre Dame, IN

#include "constitutive_model.h"
#include "construct_linearization_parts.h"
#include <ttl/ttl.h>

namespace {
  template<int R,int D = 3,class S = double>
  using Tensor = ttl::Tensor<R,D,S>;

  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;
}

/// \param[out] Kuu_a_out 
/// \param[in] A_in
/// \param[in] S_in
/// \param[in] C_in
/// \param[in] Pa_in
/// \param[in] L_in
/// \param[in] drdtaus_a
/// \return non-zero on internal error 
int compute_Kuu_a(double *Kuu_a_out, double *A_in,  double *S_in, 
                  double *C_in,      double *Pa_in, double *L_in, double drdtaus_a)
{
  Tensor<4, 3, double*> Kuu_a(Kuu_a_out);
  Tensor<2, 3, double*> A(A_in);
  Tensor<2, 3, double*> S(S_in);
  Tensor<2, 3, double*> C(C_in);
  Tensor<2, 3, double*> Pa(Pa_in);
  Tensor<4, 3, double*> L(L_in);
       
  int err = 0;
  Kuu_a = {};

  enum {AA,L_C,Dtau_dM,F2end};

  Tensor<2> F2[F2end];   //declare an array of tensors and initialize them to 0
  for (int a = 0; a < F2end; a++) {
    F2[a] = {};
  } 
   
  F2[AA](i,j) += Pa(i,k) * S(k,j);             //F2[AA] = 1*F2[AA] + 1*Pa*S
  F2[AA](i,j) += S(i,k) * Pa(j,k).to(k,j);     //F2[AA] = 1*F2[AA] + 1*S*Pa'
  F2[L_C](i,j) = L(i,j,k,l) * C(k,l);          //F2[L_C] = L:C
  F2[AA](i,j) += F2[L_C](i,k) * Pa(k,j);       //F2[AA] = 1*F2[AA] + 1*F2[L_C]*Pa
  F2[Dtau_dM](i,j) = A(k,i).to(i,k) * F2[AA](k,j); //F2[Dtau_dM] = A'*F2[AA]

  //Kronecker product scaled by drdtaus_a
  Kuu_a(i,j,k,l) = drdtaus_a * (F2[Dtau_dM])(i,j) * Pa(k,l);

  return err;
}

/// \param[out] Kuu_b_out
/// \param[in] MI_in
/// \return non-zero on internal error 
int compute_Kuu_b(double *Kuu_b_out, double *MI_in)
{
  int err = 0;
  
  Tensor<4, 3, double*> Kuu_b(Kuu_b_out);
  Tensor<2, 3, double*> MI(MI_in);

  Kuu_b = {};
  Kuu_b(i,j,k,l) = MI(j,i) * MI(l,k) + MI(l,i) * MI(j,k); 

 return err;
}
