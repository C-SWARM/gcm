//
//  KMS-IJSS2017-templates.hpp
//  ttl-learning
//
//  Created by Alberto Salvadori on 1/5/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

/*!
 
 \file KMS-IJSS2017.h
 \brief Template code file for the KMS-IJSS2017 model
 
 This file contains template classes and methods for the poroviscoplastic model desccribed in the paper
 "A. Krairi, K. Matouš, A. Salvadori, A poro-viscoplastic constitutive model for granular materials at finite strain, submitted to IJSS (2017)".
 
 */

#ifndef KMS_IJSS2017_implicit_templates_hpp
#define KMS_IJSS2017_implicit_templates_hpp


// include

#include <stdio.h>
#include <sstream>
#include <random>
#include <iostream>
#include <iomanip>
#include <math.h>




// temporary operations on Tensors
namespace{
  
#include "ttl-tools.h"
  
}


// class KMS_IJSS2017_Implicit_BE
// backward euler integrator
// ******************************


template <int dim>
KMS_IJSS2017_Implicit_BE<dim>::KMS_IJSS2017_Implicit_BE(KMS_IJSS2017_Parameters* P, KMS_IJSS2017_IO* InpOut, TimeIntegrationDataManager* TIDM )
: KMS_IJSS2017_Integration_Algorithms<dim>( P, InpOut, TIDM )
//! Constructor with material parameters and IO
{}


// IO - methods
// ------------

template <int dim>
void KMS_IJSS2017_Implicit_BE<dim>::AsAString( std::string& str, bool Verbose )
//! This method prints the model features as a string
{
  
  if ( Verbose )
  {
    // Parameters
    if (this->Parameters == NULL )
      str += " No parameters have been set, yet. \n";
    else
      this->Parameters->AsAString(str);
    
    // Integrators
    str += "\n";
    str += "KMS_IJSS2017 Implicit Integrator for dimension ";
    str += std::to_string(dim);
    str += ". \n";
  }
  
}


// Protected Methods
// -----------------

template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DnormdevtauDtau( const FTensors& A, const bool ZeroWarning )
//! Derivative of the Frobenus norm of dev(A) wrt to A
{
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  
  double den = FrobeniusNorm( dev( A ) );
  
  if ( std::fpclassify(den) == FP_ZERO )
  {
    if ( ZeroWarning )
    {
      std::cerr << " WARNING: KMS_IJSS2017_Implicit_BE<dim>::DtaubarDtau( const FTensors& A, const bool ZeroWarning ) ";
      std::cerr << " Derivative of the deviatoric part not defined since || dev( A ) || = " << den << ". \n";
      std::cerr << " Code did not abort, returned ZERO but outcomes might be wrong. \n";
    }
    
    return A(i,j) *  0;
  }
  
  double oneoverden = 1.0 / den;
  
  ttl::Tensor<2, dim, double> devA;
  devA(i,j) = dev( A )(i,j) * oneoverden;
  
  return devA ;
  
}



template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DSDM( const FTensors& F, const FTensors& eF, const FTensors& pF, const double pc )
//! Derivative of the Tensor S (Second Piola-Kirchoff stress) wrt to M.
//! This function returns the derivative of S wrt to M, at a given F, Fe, Fp.
//! Note that these three tensors are not related in the integration algorithm since F is known at step n+1, whereas
//! Fe and Fp are known only at step n and iteration k and so does pc.
{
  
  // indexes
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'e'> e;
  static constexpr ttl::Index<'f'> f;
  static constexpr ttl::Index<'g'> g;
  static constexpr ttl::Index<'h'> h;
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'z'> z;
  
  FTensors pFinv;
  pFinv = ttl::inverse( pF );
  
  double Je = det( eF ), onethird = 1.0/3.0, powJe =  pow(Je, - 2.0 * onethird ) ;
  double coh = this->c( pc );
  double alpha = pow( 1.0 -  coh / this->Parameters->c_inf, this->Parameters->pl_n ) / this->Parameters->K_kappa  ;

  // delta
  ttl::Tensor<4, dim, double> Id4 = ttl::identity( i,j,k,l ); // ttlIdentity<dim,double>(); // ttl::Delta<4,dim,double>(1); //
  ttl::Tensor<2, dim, double> Id2 = ttl::identity( i,j ); // ttl::Delta<2,dim,double>(1);
  
  // Left Caucgy-Green tensor ...
  FTensors Ce, InvCe;
  Ce(e,f) = eF(z,e) * eF(z,f);
  InvCe = ttl::inverse( Ce );
  
  // ... and its derivative
  ttl::Tensor<4, dim, double> DCeDeF;
  DCeDeF(e,f,g,h) = Id4(z,e,g,h) * eF(z,f) + eF(z,e) * Id4(z,f,g,h);
  
  // ... another derivative
  ttl::Tensor<4, dim, double> DeFDpFInv;
  DeFDpFInv(g,h,i,j) = F(g,z) * Id4(z,h,i,j);
  
  // ... and yet another
  ttl::Tensor<4, dim, double> DpFInvDM;
  DpFInvDM(i,j,l,k) = pFinv(i,z) * Id4(z,j,l,k);
  
  // derivatives of the Second Piola-Kirchoff stress
  ttl::Tensor<4, dim, double> DSisoDCe;
  DSisoDCe(a,b,e,f) =
   this->shearmodulus( pc ) * ( - onethird * powJe *  InvCe(e,f) ) * ( Id2(a,b) - Trace( Ce ) * onethird * InvCe(a,b) ) +
   this->shearmodulus( pc ) * powJe * onethird * ( -  Id2(e,f) * InvCe(a,b) + Trace( Ce ) * InvCe(a,e) * InvCe(f,b) );

  ttl::Tensor<4, dim, double> DSvolDCe;
  DSvolDCe(a,b,e,f) =
   ( this->bulkmodulus( pc ) + alpha * pow( Je, -alpha-1 ) * ( coh + this->Parameters->K_p0 ) + 0.5 * this->bulkmodulus( pc ) * log( Je ) ) * 0.5 * Je * InvCe(e,f) * InvCe(a,b)
  -( coh + 0.5 *  this->bulkmodulus( pc ) * ( Je -1 ) - pow( Je, -alpha ) * ( coh + this->Parameters->K_p0 ) + 0.5 * this->bulkmodulus( pc ) * Je * log( Je )  ) * InvCe(a,e) * InvCe(f,b);
  
  ttl::Tensor<4, dim, double> DSDCe;
  DSDCe(a,b,e,f) = DSisoDCe(a,b,e,f) + DSvolDCe(a,b,e,f);
  
  return DSDCe(a,b,e,f) * DCeDeF(e,f,g,h) * DeFDpFInv(g,h,i,j) * DpFInvDM(i,j,l,k);

}




template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DSisoDM( const FTensors& F, const FTensors& eF, const FTensors& pF, const double pc )
//! Derivative of the Tensor S (Second Piola-Kirchoff stress) wrt to M.
//! This function returns the derivative of S wrt to M, at a given F, Fe, Fp.
//! Note that these three tensors are not related in the integration algorithm since F is known at step n+1, whereas
//! Fe and Fp are known only at step n and iteration k and so does pc.
{
  
  // debug
  std::cout << "In DSisoDM. \n";
  std::string outstr="", localstr = "";

  // indexes
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'e'> e;
  static constexpr ttl::Index<'f'> f;
  static constexpr ttl::Index<'g'> g;
  static constexpr ttl::Index<'h'> h;
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'x'> x;
  
  FTensors pFinv;
  pFinv = ttl::inverse( pF );

  PrintInVectorForm( pFinv, outstr );
  std::cout << " pFinv = " << outstr ;
  outstr="";
  
  double Je = det( eF ), onethird = 1.0/3.0, powJe =  pow(Je, - 2.0 * onethird ) ;
  double coh = this->c( pc );
  double alpha = pow( 1.0 -  coh / this->Parameters->c_inf, this->Parameters->pl_n ) / this->Parameters->K_kappa  ;
  
  std::cout << " Je = " << Je << " \n";
  std::cout << " coh = " << coh << " \n";
  std::cout << " alpha = " << alpha << " \n";

  // delta
  ttl::Tensor<4, dim, double> Id4 = ttl::identity( i,j,k,l ); // ttlIdentity<dim,double>(); // ttl::Delta<4,dim,double>(1); //
  ttl::Tensor<2, dim, double> Id2 = ttl::identity( i,j ); // ttl::Delta<2,dim,double>(1);
  
  // Right Cauchy-Green tensor ...
  FTensors Ce, InvCe, Temp;
  Ce(e,f) = eF(z,e) * eF(z,f);
  InvCe = ttl::inverse( Ce );
  
  Temp(i,j) = ttl::identity(i,j,k,l) * Ce(k,l);
  PrintInVectorForm( Ce, outstr );
  std::cout << " Ce = " << outstr;
  outstr="";
  PrintInVectorForm( Temp, outstr );
  std::cout << " Id4(i,j,k,l) * Ce(k,l) = " << outstr;
  outstr="";
  
  PrintInVectorForm( InvCe, outstr );
  std::cout << " InvCe = " << outstr;
  outstr="";

  PrintInVectorForm( eF, outstr );
  std::cout << " eF = " << outstr;
  outstr="";

  PrintInVectorForm( F, outstr );
  std::cout << " F = " << outstr;
  outstr="";

  // ... and its derivative
  ttl::Tensor<4, dim, double> DCeDeF;
  DCeDeF(e,f,g,h) = Id4(z,e,g,h) * eF(z,f) + eF(x,e) * Id4(x,f,g,h);
  
  // scrivere DCeDeF che e` un tensore del quarto ordine
  for ( unsigned i=0; i<dim; i++ )
    for ( unsigned j=0; j<dim; j++ )
    {
      localstr = "";
      PrintInVectorForm(DCeDeF, i, j, localstr);
      std::cout << ". i =" << i << ", j =" << j << ", DCeDeF(i,j) = " << localstr;
    }
  std::cout << "\n";
  
  // ... another derivative
  ttl::Tensor<4, dim, double> DeFDpFInv;
  DeFDpFInv(g,h,i,j) = F(g,z) * Id4(z,h,i,j);
  
  // scrivere DeFDpFInv che e` un tensore del quarto ordine
  for ( unsigned i=0; i<dim; i++ )
    for ( unsigned j=0; j<dim; j++ )
    {
      localstr = "";
      PrintInVectorForm(DeFDpFInv, i, j, localstr);
      std::cout << ". i =" << i << ", j =" << j << ", DeFDpFInv(i,j) = " << localstr;
    }
  std::cout << "\n";

  // ... and yet another
  ttl::Tensor<4, dim, double> DpFInvDM;
  DpFInvDM(i,j,l,k) = pFinv(i,z) * Id4(z,j,l,k);
  
  // scrivere DpFInvDM che e` un tensore del quarto ordine
  for ( unsigned i=0; i<dim; i++ )
    for ( unsigned j=0; j<dim; j++ )
    {
      localstr = "";
      PrintInVectorForm(DpFInvDM, i, j, localstr);
      std::cout << ". i =" << i << ", j =" << j << ", DpFInvDM(i,j) = " << localstr;
    }
  std::cout << "\n";

  // derivatives of the Second Piola-Kirchoff stress
  ttl::Tensor<4, dim, double> DSisoDCe;
  DSisoDCe(a,b,e,f) =
  this->shearmodulus( pc ) * ( - onethird * powJe *  InvCe(e,f) ) * ( Id2(a,b) - Trace( Ce ) * onethird * InvCe(a,b) ) +
  this->shearmodulus( pc ) * powJe * onethird * ( - Id2(e,f) * InvCe(a,b) + Trace( Ce ) * InvCe(a,e) * InvCe(f,b) );
  
  // scrivere DSisoDCe che e` un tensore del quarto ordine
  for ( unsigned i=0; i<dim; i++ )
    for ( unsigned j=0; j<dim; j++ )
    {
      localstr = "";
      PrintInVectorForm(DSisoDCe, i, j, localstr);
      std::cout << ". i =" << i << ", j =" << j << ", DSisoDCe(i,j) = " << localstr;
    }
  std::cout << "\n";

  ttl::Tensor<4, dim, double> DSDCe;
  DSDCe(a,b,e,f) = DSisoDCe(a,b,e,f);
  
  
  std::cout << "In DSisoDM, end. \n\n";
  
  
  return DSDCe(a,b,e,f) * DCeDeF(e,f,g,h) * DeFDpFInv(g,h,i,j) * DpFInvDM(i,j,l,k);
  
}


template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DSvolDM( const FTensors& F, const FTensors& eF, const FTensors& pF, const double pc )
//! Derivative of the Tensor S (Second Piola-Kirchoff stress) wrt to M.
//! This function returns the derivative of S wrt to M, at a given F, Fe, Fp.
//! Note that these three tensors are not related in the integration algorithm since F is known at step n+1, whereas
//! Fe and Fp are known only at step n and iteration k and so does pc.
{
  
  // indexes
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'e'> e;
  static constexpr ttl::Index<'f'> f;
  static constexpr ttl::Index<'g'> g;
  static constexpr ttl::Index<'h'> h;
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'z'> z;
  
  FTensors pFinv;
  pFinv = ttl::inverse( pF );
  
  double Je = det( eF ); // , onethird = 1.0/3.0, powJe =  pow(Je, - 2.0 * onethird ) ;
  double coh = this->c( pc );
  double alpha = pow( 1.0 -  coh / this->Parameters->c_inf, this->Parameters->pl_n ) / this->Parameters->K_kappa  ;
  
  // delta
  ttl::Tensor<4, dim, double> Id4 = ttl::identity( i,j,k,l ); // ttlIdentity<dim,double>(); // ttl::Delta<4,dim,double>(1); //
  ttl::Tensor<2, dim, double> Id2 = ttl::identity( i,j ); // ttl::Delta<2,dim,double>(1);
  
  // Left Caucgy-Green tensor ...
  FTensors Ce, InvCe;
  Ce(e,f) = eF(z,e) * eF(z,f);
  InvCe = ttl::inverse( Ce );
  
  // ... and its derivative
  ttl::Tensor<4, dim, double> DCeDeF;
  DCeDeF(e,f,g,h) = Id4(z,e,g,h) * eF(z,f) + eF(z,e) * Id4(z,f,g,h);
  
  // ... another derivative
  ttl::Tensor<4, dim, double> DeFDpFInv;
  DeFDpFInv(g,h,i,j) = F(g,z) * Id4(z,h,i,j);
  
  // ... and yet another
  ttl::Tensor<4, dim, double> DpFInvDM;
  DpFInvDM(i,j,l,k) = pFinv(i,z) * Id4(z,j,l,k);
  
  ttl::Tensor<4, dim, double> DSvolDCe;
  DSvolDCe(a,b,e,f) =
  ( this->bulkmodulus( pc ) + alpha * pow( Je, -alpha-1 ) * ( coh + this->Parameters->K_p0 ) + 0.5 * this->bulkmodulus( pc ) * log( Je ) ) * 0.5 * Je * InvCe(e,f) * InvCe(a,b)
  -( coh + 0.5 *  this->bulkmodulus( pc ) * ( Je -1 ) - pow( Je, -alpha ) * ( coh + this->Parameters->K_p0 ) + 0.5 * this->bulkmodulus( pc ) * Je * log( Je )  )  * 0.5 * ( InvCe(a,e) * InvCe(b,f) + InvCe(b,e) * InvCe(a,f) ); // * InvCe(a,e) * InvCe(f,b); //
  
  ttl::Tensor<4, dim, double> DSDCe;
  DSDCe(a,b,e,f) = DSvolDCe(a,b,e,f);
  
  // debug
  std::cout << "In DSisoDM. \n";
  std::cout << " Je = " << Je << " \n";
  std::cout << " coh = " << coh << " \n";
  std::cout << " alpha = " << alpha << " \n";
  
  std::cout << "In DSisoDM, end. \n\n";

  return DSDCe(a,b,e,f) * DCeDeF(e,f,g,h) * DeFDpFInv(g,h,i,j) * DpFInvDM(i,j,l,k);
  
}





template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DGammaDM( const ttl::Tensor<4, dim, double>& dsdm, const FTensors& M, const FTensors& S)
//! Derivative of the Tensor Gamma wrt to M
//! Tensor Gamma is defined as ( det M )^(-1) * M * S(M) * M
//! This function returns the derivative of Gamma wrt to M, at a given M and S
//! Note that DSDM is passed as an argument, i.e. it has to be evaluated beforehand
{
  
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;

  ttl::Tensor<4, dim, double> Der;
  ttl::Tensor<4, dim, double> Id4 = ttl::identity( i,j,k,l ); // ttlIdentity<dim,double>(); // ttl::Delta<4,dim,double>(1); //
  
  // In matrix calculus, Jacobi's formula expresses the derivative of the determinant of a matrix A in terms of the adjugate of A.
  // d Det[A] / dA_ij = adjT(A)_ij
  // the tensor adjTM is the transpose of the adjugate of M
  ttl::Tensor<2, dim, double> adjTM =  adjugateTranspose( M );
 
  // the derivative of ( det M )^(-1) * M * S(M) * M has 4 contributions
  // listed below. See also the accompanying notes.
  // 1
  double detMm1 = 1.0 / det(M);
  Der(z,w,l,k) = - adjTM(l,k) * M(z,a) * S(a,b) * M(w,b) * detMm1 * detMm1;
  // 2
  Der(z,w,l,k) += detMm1 * Id4(z,a,l,k) * S(a,b) * M(w,b);
  // 3
  Der(z,w,l,k) += detMm1 * M(z,a) * dsdm(a,b,l,k) * M(w,b);
  // 4
  Der(z,w,l,k) += detMm1 * M(z,a) * S(a,b) * Id4(w,b,l,k);
  
  return Der;
  
}



template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::HattedDGammaDM( const ttl::Tensor<4, dim, double>& dsdm, const double dHpdpc, const FTensors& M, const FTensors& S)
//! Derivative of the Tensor Gamma wrt to M
//! Tensor Gamma is defined as ( det M )^(-1) * M * S(M) * M
//! This function returns the derivative of Gamma wrt to M, at a given M and S
//! Note that DSDM is passed as an argument, i.e. it has to be evaluated beforehand
//!
//! See page 2, "debugging and code optimization" notes
{
  
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  
  ttl::Tensor<4, dim, double> Der = {};
  ttl::Tensor<4, dim, double> Id4 = ttl::identity( i,j,k,l ); // ttlIdentity<dim,double>(); // ttl::Delta<4,dim,double>(1); //
  
  // In matrix calculus, Jacobi's formula expresses the derivative of the determinant of a matrix A in terms of the adjugate of A.
  // d Det[A] / dA_ij = adjT(A)_ij
  // the tensor adjTM is the transpose of the adjugate of M
  ttl::Tensor<2, dim, double> adjTM =  adjugateTranspose( M );
  
  // the derivative of ( det M )^(-1) * M * S(M) * M has 4 contributions
  // listed below. See also the accompanying notes.
  // 1
  double detMm1 = 1.0 / det(M);
  Der(z,w,l,k) = - dHpdpc * adjTM(l,k) * M(z,a) * S(a,b) * M(w,b) * detMm1 * detMm1;
  // 2
  Der(z,w,l,k) += dHpdpc * detMm1 * Id4(z,a,l,k) * S(a,b) * M(w,b);
  // 3 - dsdm already contains dHpdpc manipulation,
  // See page 2, "debugging and code optimization" notes
  Der(z,w,l,k) += detMm1 * M(z,a) * dsdm(a,b,l,k) * M(w,b);
  // 4
  Der(z,w,l,k) += dHpdpc * detMm1 * M(z,a) * S(a,b) * Id4(w,b,l,k);
  
  return Der;
  
}



template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DKSDM( const ttl::Tensor<4, dim, double>& dgammadm, const FTensors& F, const FTensors& pF)
//! Derivative of the Kirchoff stress wrt to M.
//! This function returns the derivative of KS wrt to M, at a given F, Fp.
//! Note that these two tensors are not related in the integration algorithm since F is known at step n+1, whereas
//! Fp is known only at step n and iteration k and so does pc.
{
  
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'x'> x;
  
  FTensors Aa, InvFp = ttl::inverse( pF );
  Aa(a,z) =  F(a,x) * InvFp(x,z);
  
  ttl::Tensor<4, dim, double> Der;
  Der(a,b,l,k) = det(pF) * Aa( a,z ) * dgammadm( z,w,l,k ) * Aa( b, w );
  
  return Der;

}


template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DgammadotdDM(const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& Mr, const double pcr )
//! Derivative of gammadotd wrt to M - see notes
//! pFn is the plastic tensor at step n. Not to be confused with the iteration r of the NR scheme. pFn is introduced in the definition of M and shall not be updated
//!     until the convergence of the NR is achieved and step n updated by step n+1 ( by method StepClose() )
//! Fnp1 is the total strain tensor at step n+1. it remains unchanged with the NR iterations on r
//! eFr is the elastic strain tensor at ITERATION r. It eventually updates the Second Piola-Kirchoff stress.
//! Sr is the Second Piola-Kirchoff stress tensor at ITERATION r.
//! Mr is equal to ( pFn . Inv(pFnp1)r ) at ITERATION r.
//! pcr is the internal variable at ITERATION r.

{
  
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  
  
  // Jp = Log(Fp) at iteration r
  double pJr = det( pFn ) / det ( Mr );
  
  // Kirchoff stress at iteration r
  FTensors KSr;
  KSr( a,b ) = pJr * eFr( a,x ) * Sr( x,z ) * eFr( b,z );
  
  //! Derivative of the Frobenus norm of dev( Kirchoff stress ) wrt to Kirchoff stress
  FTensors dtaubardtau = DnormdevtauDtau( KSr );
  
  //! Derivative of the Kirchoff stress wrt to M.
  ttl::Tensor<4, dim, double> dsdm = DSDM( Fnp1, eFr, pFn, pcr );
  ttl::Tensor<4, dim, double> dgammadm = DGammaDM( dsdm, Mr, Sr );
  ttl::Tensor<4, dim, double> dtaudm = DKSDM( dgammadm, Fnp1, pFn );
  
  //! Derivative of the Frobenius norm of dev( Kirchoff stress ) wrt to M
  FTensors dtaubardM, temp;
  dtaubardM( l,k ) = dtaubardtau( a,b ) * dtaudm ( a,b,l,k );
 
  // taubar is the Frobenius norm of its deviatoric part
  double taubar = FrobeniusNorm( dev( KSr ) );
  
  // model parameters at iteration r
  double dr = this->d( pcr );
  double gammataur = this->g_tau( pcr );
  double oneoverm = 1.0 / this->Parameters->flr_m ;
  
  double coeff = this->Parameters->flr_gamma0 * ( 1.0 - 1.0 / dr ) * oneoverm  / gammataur  *  pow( taubar / gammataur , oneoverm - 1 ) ;
  
  temp( l,k ) = coeff * dtaubardM( l,k ) ;
  
  return temp;
}


template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DgammadotdDM( const ttl::Tensor<4, dim, double>& dtaudm, const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& KSr, const FTensors& Mr, const double pcr )
//! Derivative of gammadotd wrt to M - see notes
//! pFn is the plastic tensor at step n. Not to be confused with the iteration r of the NR scheme. pFn is introduced in the definition of M and shall not be updated
//!     until the convergence of the NR is achieved and step n updated by step n+1 ( by method StepClose() )
//! Fnp1 is the total strain tensor at step n+1. it remains unchanged with the NR iterations on r
//! eFr is the elastic strain tensor at ITERATION r. It eventually updates the Second Piola-Kirchoff stress.
//! Sr is the Second Piola-Kirchoff stress tensor at ITERATION r.
//! Mr is equal to ( pFn . Inv(pFnp1)r ) at ITERATION r.
//! pcr is the internal variable at ITERATION r.

{
  
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  
  //! Derivative of the Frobenus norm of dev( Kirchoff stress ) wrt to Kirchoff stress
  FTensors dtaubardtau = DnormdevtauDtau( KSr );
  
  //! Derivative of the Frobenius norm of dev( Kirchoff stress ) wrt to M
  FTensors dtaubardM, temp;
  dtaubardM( l,k ) = dtaubardtau( a,b ) * dtaudm ( a,b,l,k );
  
  // taubar is the Frobenius norm of its deviatoric part
  double taubar = FrobeniusNorm( dev( KSr ) );
  
  // model parameters at iteration r
  double dr = this->d( pcr );
  double gammataur = this->g_tau( pcr );
  double oneoverm = 1.0 / this->Parameters->flr_m ;
  
  double coeff = this->Parameters->flr_gamma0 * ( 1.0 - 1.0 / dr ) * oneoverm  / gammataur  *  pow( taubar / gammataur , oneoverm - 1 ) ;
  
  temp( l,k ) = coeff * dtaubardM( l,k ) ;
  
  return temp;
}


template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DgammadotvDM(const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& KSr, const FTensors& Mr, const double pcr )
//! Derivative of gammadotd wrt to M - see notes
//! pFn is the plastic tensor at step n. Not to be confused with the iteration r of the NR scheme. pFn is introduced in the definition of M and shall not be updated
//!     until the convergence of the NR is achieved and step n updated by step n+1 ( by method StepClose() )
//! Fnp1 is the total strain tensor at step n+1. it remains unchanged with the NR iterations on r
//! eFr is the elastic strain tensor at ITERATION r. It eventually updates the Second Piola-Kirchoff stress.
//! Sr is the Second Piola-Kirchoff stress tensor at ITERATION r.
//! KSr is the Kirchoff stress tensor at ITERATION r.
//! Mr is equal to ( pFn . Inv(pFnp1)r ) at ITERATION r.
//! pcr is the internal variable at ITERATION r.

{
  
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  
  
  
  // Kirchoff pressure. Beware of the sign.
  double pi = - Trace( KSr ) / 3.0;
  
  // model parameters at iteration r
  double pim = this->pi_m( pcr );
  double gpir = this->g_pi( pi, pcr );
  double oneoverm = 1.0 / this->Parameters->flr_m ;
  
  // if pi is not larger than pim return zero. This check should be done before invoking
  // this function, though...
  if ( pi <= pim )
    return {};
  
  //! Derivative of the Kirchoff pressure wrt to Kirchoff stress. Beware of the sign.
  // ttl::Tensor<2, dim, double> Id2 = ttl::identity( a,b ); // ttl::Delta<2,dim,double>();
  FTensors dpidtau =  ( - 1.0 / 3.0 ) * ttl::identity<2>( a,b ); //* Id2( a,b );
  
  //! Derivative of the Kirchoff stress wrt to M.
  ttl::Tensor<4, dim, double> dsdm = DSDM( Fnp1, eFr, pFn, pcr );
  ttl::Tensor<4, dim, double> dgammadm = DGammaDM( dsdm, Mr, Sr );
  ttl::Tensor<4, dim, double> dtaudm = DKSDM( dgammadm, Fnp1, pFn );
  
  //! Derivative of the Frobenius norm of dev( Kirchoff stress ) wrt to M
  FTensors dpidM, temp;
  dpidM( l,k ) = dpidtau( a,b ) * dtaudm ( a,b,l,k );
  
  double coeff = this->Parameters->flr_gamma0 * oneoverm  / gpir  *  pow( ( pi - pim ) / gpir , oneoverm - 1 ) ;
  
  temp( l,k ) = coeff * dpidM( l,k ) ;
  
  return temp;
}




template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DgammadotvDM( const ttl::Tensor<4, dim, double>& dtaudm, const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& KSr, const FTensors& Mr, const double pcr )
//! Derivative of gammadotd wrt to M - see notes
//! pFn is the plastic tensor at step n. Not to be confused with the iteration r of the NR scheme. pFn is introduced in the definition of M and shall not be updated
//!     until the convergence of the NR is achieved and step n updated by step n+1 ( by method StepClose() )
//! Fnp1 is the total strain tensor at step n+1. it remains unchanged with the NR iterations on r
//! eFr is the elastic strain tensor at ITERATION r. It eventually updates the Second Piola-Kirchoff stress.
//! Sr is the Second Piola-Kirchoff stress tensor at ITERATION r.
//! KSr is the Kirchoff stress tensor at ITERATION r.
//! Mr is equal to ( pFn . Inv(pFnp1)r ) at ITERATION r.
//! pcr is the internal variable at ITERATION r.

{
  
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  
  
  // Kirchoff pressure. Beware of the sign.
  double pi = - Trace( KSr ) / 3.0;
  
  // model parameters at iteration r
  double pim = this->pi_m( pcr );
  double gpir = this->g_pi( pi, pcr );
  double oneoverm = 1.0 / this->Parameters->flr_m ;

  // if pi is not larger than pim return zero. This check should be done before invoking
  // this function, though...
  if ( pi <= pim )
    return {};
  
  //! Derivative of the Kirchoff pressure wrt to Kirchoff stress. Beware of the sign.
  ttl::Tensor<2, dim, double> Id2 = ttl::identity( a,b ); // ttl::Delta<2,dim,double>();
  FTensors dpidtau = Id2( a,b ) * ( - 1.0 / 3.0 );
  
  //! Derivative of the Frobenius norm of dev( Kirchoff stress ) wrt to M
  FTensors dpidM, temp;
  dpidM( l,k ) = dpidtau( a,b ) * dtaudm ( a,b,l,k );
  
  double coeff = this->Parameters->flr_gamma0 * oneoverm  / gpir  *  pow( ( pi - pim ) / gpir , oneoverm - 1 ) ;
  
  temp( l,k ) = coeff * dpidM( l,k ) ;
  
  return temp;
}






template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DPsidDM( const FTensors& devS, const FTensors& F, const FTensors& eF, const FTensors& pF, const double pc )
//! Derivative of the Tensor S (Second Piola-Kirchoff stress) wrt to M.
//! This function returns the derivative of S wrt to M, at a given F, Fe, Fp.
//! Note that these three tensors are not related in the integration algorithm since F is known at step n+1, whereas
//! Fp and Fe are known only at step n and iteration r and so does pc.
//! devS is the deviator of the Second Piola-Kirchoff stress, at iteration r
{
  
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'c'> c;
  static constexpr ttl::Index<'d'> d;
  static constexpr ttl::Index<'e'> e;
  static constexpr ttl::Index<'f'> f;

  // tensor definitions and initialization
  ttl::Tensor<4, dim, double> dpsiddm {};
  ttl::Tensor<2, dim, double> Id2 = ttl::identity( i,c ); // ttl::Delta<2,dim,double>();  // ttlIdentity<dim,double>();
  
  // norm of the deviatoric part of the Second-Piola Kirchoff stress
  double normPsid = FrobeniusNorm( devS );
 
  if ( std::fpclassify( normPsid ) == FP_ZERO )
  {
    // this test should not be necessary, since DPsidDM shall not be invoked if devS = 0
    // in case though, a zero-tensor is returned
  }
  else
  {
    double den = pow( normPsid, -3.0 );
    
    dpsiddm( i,j,l,k ) = den
     * ( - devS( i,j ) * devS( c,d ) + normPsid * normPsid * Id2( i,c ) * Id2( j,d ) )
     * ( Id2( c,e ) * Id2( d,f ) - 1.0 / 3.0 * Id2( c,d ) * Id2( e,f ) )
     * DSDM( F, eF, pF, pc )( e,f,l,k );
  }
  
  return dpsiddm;
  
}
  
  
  
template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DRMDM( const double Deltat, const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& Mr, const double pcr )
//! Diagonal contribution of the residual wrt to M
//! this implementation is supposed to be readable, although not optimized for speed
//!
//! pFn is the plastic tensor at step n. Not to be confused with the iteration r of the NR scheme. pFn is introduced in the definition of M and shall not be updated
//!     until the convergence of the NR is achieved and step n updated by step n+1 ( by method StepClose() )
//! Fnp1 is the total strain tensor at step n+1. it remains unchanged with the NR iterations on r
//! eFr is the elastic strain tensor at ITERATION r. It eventually updates the Second Piola-Kirchoff stress.
//! Sr is the Second Piola-Kirchoff stress tensor at ITERATION r.
//! Mr is equal to ( pFn . Inv(pFnp1)r ) at ITERATION r.
//! pcr is the internal variable at ITERATION r.
//!
//! dRMdM has the following four contributions:
//! 1 - D gammadotd / DM \times Psid
//! 2 - gammadotd \times  D Psid / DM
//! 3 - D gammadotv / DM \times Psiv
//! 4 - D M / DM - identity operator
{
  
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  
  // definition of drdm with initialization to zero
  ttl::Tensor<4, dim, double> drdm = {};
 
  // evaluation of the fourth order tensor dtaudm
  
  // Jp = Log(Fp) at iteration r
  double pJr = det( pFn ) / det ( Mr );
  
  // Kirchoff stress at iteration r
  FTensors KSr;
  KSr( a,b ) = pJr * eFr( a,x ) * Sr( x,z ) * eFr( b,z );
    
  //! Derivative of the Kirchoff stress wrt to M.
  ttl::Tensor<4, dim, double> dsdm = DSDM( Fnp1, eFr, pFn, pcr );
  ttl::Tensor<4, dim, double> dgammadm = DGammaDM( dsdm, Mr, Sr );
  ttl::Tensor<4, dim, double> dtaudm = DKSDM( dgammadm, Fnp1, pFn );
  
  // 1 and 2
  
  // deviatoric part of the Second-Piola Kirchoff stress
  FTensors Psid = dev( Sr );
  
  // and its norm
  double normPsid = FrobeniusNorm( Psid );

  // if Psid is zero, there is no sense to estimate dgammadotddm since Psid contracted
  // to DeltaM will always provide ZERO. Therefore, in such a case the contribution of
  // D gammadotd / DM \times Psid to DR/DM is just ZERO
  
  if ( std::fpclassify( normPsid ) == FP_ZERO )
  { }
  else
  {
    
    // norm dev( KSr )
    double taubar =  FrobeniusNorm( dev( KSr ) );
    
    // DPsidDM is estimated with Psid = dev( Sr )
    ttl::Tensor<4, dim, double> dpsiddm = DPsidDM( Psid, Fnp1, eFr, pFn, pcr );
    
    // DgammadotdDM
    FTensors dgammadotddm = DgammadotdDM( dtaudm, pFn, Fnp1, eFr, Sr, KSr, Mr, pcr );
    
    // diadic product contribution
    drdm( i,j,l,k ) = Deltat * ( Psid( i,j ) * ( 1.0 / normPsid ) * dgammadotddm( l,k ) + this->gammadot_d( taubar, pcr ) * dpsiddm( i,j,l,k ) );
    
  }
  
  // 3
  
  // Kirchoff pressure. Beware of the sign.
  double pi = - Trace( KSr ) / 3.0;
  double pim = this->pi_m( pcr );
  
  if (  ( pi < pim ) || ( pcr >= this->Parameters->cf_pcinf ) )
  { }
  else
  {
    // volumetric part of flow rule tensor
    
    ttl::Tensor<2, dim, double> Id2 = ttl::identity( i,j ); // ttl::Delta<2,dim,double>();
    FTensors Psiv = Id2( i,j ) * ( - this->betaC( pcr ) / 3.0 );
    
    // DgammadotvdDM(const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& KSr, const FTensors& Mr, const double pcr )
    FTensors dgammadotvdm = DgammadotvDM( dtaudm, pFn, Fnp1, eFr, Sr, KSr, Mr, pcr );
    
    // diadic product contribution
    drdm( i,j,l,k ) = drdm( i,j,l,k ) + Deltat * ( Psiv( i,j ) * dgammadotvdm( l,k ) );
  }
  
  // 4
  
  // DMDM
  ttl::Tensor<4, dim, double> dmdm = ttl::identity( i,j,l,k );
  
  // sum
  drdm( i,j,l,k ) = drdm( i,j,l,k ) + dmdm( i,j,l,k );
  
  return drdm;
  
  
}



template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::VerboseDRMDM( const double Deltat, const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& Mr, const double pcr )
//! Diagonal contribution of the residual wrt to M
//! this implementation is supposed to be readable, although not optimized for speed
//!
//! pFn is the plastic tensor at step n. Not to be confused with the iteration r of the NR scheme. pFn is introduced in the definition of M and shall not be updated
//!     until the convergence of the NR is achieved and step n updated by step n+1 ( by method StepClose() )
//! Fnp1 is the total strain tensor at step n+1. it remains unchanged with the NR iterations on r
//! eFr is the elastic strain tensor at ITERATION r. It eventually updates the Second Piola-Kirchoff stress.
//! Sr is the Second Piola-Kirchoff stress tensor at ITERATION r.
//! Mr is equal to ( pFn . Inv(pFnp1)r ) at ITERATION r.
//! pcr is the internal variable at ITERATION r.
//!
//! dRMdM has the following four contributions:
//! 1 - D gammadotd / DM \times Psid
//! 2 - gammadotd \times  D Psid / DM
//! 3 - D gammadotv / DM \times Psiv
//! 4 - D M / DM - identity operator
{
  
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  
  
  // definition of drdm with initialization to zero
  ttl::Tensor<4, dim, double> drdm = {};

  // 1 and 2
  
  // deviatoric part of the Second-Piola Kirchoff stress
  FTensors Psid = dev( Sr );
  
  // and its norm
  double normPsid = FrobeniusNorm( Psid );
  
  // if Psid is zero, there is no sense to estimate dgammadotddm since Psid contracted
  // to DeltaM will always provide ZERO. Therefore, in such a case the contribution of
  // D gammadotd / DM \times Psid to DR/DM is just ZERO
  
  if ( std::fpclassify( normPsid ) == FP_ZERO )
  { }
  else
  {
    
    // Kirchoff stress at iteration r
    FTensors KSr;
    KSr( a,b ) = det( pFn ) / det( Mr ) * eFr( a,x ) * Sr( x,z ) * eFr( b,z );
    
    // norm dev( KSr )
    double taubar =  FrobeniusNorm( dev( KSr ) );
    
    // DPsidDM is estimated with Psid = dev( Sr )
    ttl::Tensor<4, dim, double> dpsiddm = DPsidDM( Psid, Fnp1, eFr, pFn, pcr );

    // DgammadotdDM
    FTensors dgammadotddm = DgammadotdDM( pFn, Fnp1, eFr, Sr, Mr, pcr );
    
    // diadic product contribution
    drdm( i,j,l,k ) = Deltat * ( Psid( i,j ) * ( 1.0 / normPsid ) * dgammadotddm( l,k ) + this->gammadot_d( taubar, pcr ) * dpsiddm( i,j,l,k ) );
    
  }
  

  // 3
  
  
  // Jp = Log(Fp) at iteration r
  double pJr = det( pFn ) / det ( Mr );
  
  // Kirchoff stress at iteration r
  FTensors KSr;
  KSr( a,b ) = pJr * eFr( a,x ) * Sr( x,z ) * eFr( b,z );
  
  // Kirchoff pressure. Beware of the sign.
  double pi = - Trace( KSr ) / 3.0;
  double pim = this->pi_m( pcr );

  if (  ( pi < pim ) || ( pcr >= this->Parameters->cf_pcinf ) )
  { }
  else
  {
    // volumetric part of flow rule tensor
    
    ttl::Tensor<2, dim, double> Id2 = ttl::identity( i,j ); // ttl::Delta<2,dim,double>();
    FTensors Psiv = Id2( i,j ) * ( - this->betaC( pcr ) / 3.0 );
    
    // DgammadotvdDM(const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& KSr, const FTensors& Mr, const double pcr )
    FTensors dgammadotvdm = DgammadotvDM( pFn, Fnp1, eFr, Sr, KSr, Mr, pcr );
    
    // diadic product contribution
    drdm( i,j,l,k ) = drdm( i,j,l,k ) + Deltat * ( Psiv( i,j ) * dgammadotvdm( l,k ) );
  }

  // 4
  
  // DMDM
  ttl::Tensor<4, dim, double> dmdm = ttl::identity( i,j,l,k );
  
  // sum
  drdm( i,j,l,k ) = drdm( i,j,l,k ) + dmdm( i,j,l,k );
  
  return drdm;
}







template <int dim>
double KMS_IJSS2017_Implicit_BE<dim>::DRpcDpc( const double Deltat, const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& Mr, const double pcr )
//! Diagonal contribution of the residual wrt to M
//! this implementation is supposed to be readable, although not optimized for speed
//!
//! pFn is the plastic tensor at step n. Not to be confused with the iteration r of the NR scheme. pFn is introduced in the definition of M and shall not be updated
//!     until the convergence of the NR is achieved and step n updated by step n+1 ( by method StepClose() )
//! Fnp1 is the total strain tensor at step n+1. it remains unchanged with the NR iterations on r
//! eFr is the elastic strain tensor at ITERATION r. It eventually updates the Second Piola-Kirchoff stress.
//! Sr is the Second Piola-Kirchoff stress tensor at ITERATION r.
//! Mr is equal to ( pFn . Inv(pFnp1)r ) at ITERATION r.
//! pcr is the internal variable at ITERATION r.
//!
{
  
  
  // Material parameters
  double a1=this->Parameters->hr_a1, a2=this->Parameters->hr_a2, l1=this->Parameters->hr_Lambda1, l2=this->Parameters->hr_Lambda2;

  // Hardening function at pcr
  double Hpc = this->HardeningLaw( pcr );
  
  // DRpcDpc
  double dRpcdpc = - det( pFn ) / ( pcr * pcr * exp( Hpc ) ) * ( a1 * l1 * exp( - l1 / pcr ) + a2 * l2 * exp( - l2 / pcr ) );

  
  return dRpcdpc;

}





template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017_Implicit_BE<dim>::RM(  const double Deltat, const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& Mr, const double pcr )
//! Right hand side contribution of the residual for M
//! RM = Deltat ( gammadotd PsiD + gammadotv Psiv ) - Id + M
//!
//! pFn is the plastic tensor at step n. Not to be confused with the iteration r of the NR scheme. pFn is introduced in the definition of M and shall not be updated
//!     until the convergence of the NR is achieved and step n updated by step n+1 ( by method StepClose() )
//! Fnp1 is the total strain tensor at step n+1. it remains unchanged with the NR iterations on r
//! eFr is the elastic strain tensor at ITERATION r. It eventually updates the Second Piola-Kirchoff stress.
//! Sr is the Second Piola-Kirchoff stress tensor at ITERATION r.
//! Mr is equal to ( pFn . Inv(pFnp1)r ) at ITERATION r.
//! pcr is the internal variable at ITERATION r.
//!
{
  
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;

  FTensors Id2 = ttl::identity( i,j ); // ttl::Delta<2,dim,double>();
  FTensors buf = {};
  
  // deviatoric part of the Second-Piola Kirchoff stress
  FTensors Psid = dev( Sr );
  
  // and its norm
  double normPsid = FrobeniusNorm( Psid );
  
  // Kirchoff stress at iteration r
  FTensors KSr;
  KSr( a,b ) = det( pFn ) / det( Mr ) * eFr( a,x ) * Sr( x,z ) * eFr( b,z );


  // if Psid is zero, there is no sense to estimate dgammadotddm since Psid contracted
  // to DeltaM will always provide ZERO. Therefore, in such a case the contribution of
  // D gammadotd / DM \times Psid to DR/DM is just ZERO
  
  if ( std::fpclassify( normPsid ) == FP_ZERO )
  { }
  else
  {
    
    // norm dev( KSr )
    double taubar =  FrobeniusNorm( dev( KSr ) );
    
    // gammadotd PsiD
    double gammadotd = this->gammadot_d( taubar, pcr );
    
    buf( i,j ) = gammadotd * Psid( i,j ) / normPsid;
  }
  
  // volumetric part of flow rule tensor
  // Kirchoff pressure. Beware of the sign.
  double pi = - Trace( KSr ) / 3.0;
  double pim = this->pi_m( pcr );
  
  // if pi < pim gammadotv is just ZERO
  if (  ( pi < pim ) || ( pcr >= this->Parameters->cf_pcinf ) )
  { }
  else
  {
    
    // gammadotv Psiv

    double gammadotv = this->gammadot_v( pi, pcr );
    FTensors Psiv = ( - this->betaC( pcr ) / 3.0 ) *  Id2( i,j );

    buf( i,j ) = ( buf( i,j ) + gammadotv * Psiv( i,j ) ) * Deltat;
  }
  
  // Mr - Id2
  
  buf( i,j ) = buf( i,j ) - Id2( i,j ) + Mr( i,j );
  
  // return
  
  return buf;
  
}



template <int dim>
double KMS_IJSS2017_Implicit_BE<dim>::Rpc(  const double Deltat, const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& Mr, const double pcr )
//! Right hand side contribution of the residual for M
//! Rpc = det( M ) - det( pFn) g(pc)
//!
//! pFn is the plastic tensor at step n. Not to be confused with the iteration r of the NR scheme. pFn is introduced in the definition of M and shall not be updated
//!     until the convergence of the NR is achieved and step n updated by step n+1 ( by method StepClose() )
//! Fnp1 is the total strain tensor at step n+1. it remains unchanged with the NR iterations on r
//! eFr is the elastic strain tensor at ITERATION r. It eventually updates the Second Piola-Kirchoff stress.
//! Sr is the Second Piola-Kirchoff stress tensor at ITERATION r.
//! Mr is equal to ( pFn . Inv(pFnp1)r ) at ITERATION r.
//! pcr is the internal variable at ITERATION r.
//!
{

  return det(Mr) - det(pFn) / exp( this->HardeningLaw( pcr ) );
  
}




// Testing Methods
// ---------------



template <int dim>
void KMS_IJSS2017_Implicit_BE<dim>::TestSuite( std::string& str, const FTensors& F0, const FTensors& F1 )
//! This method tests some of the features of the KMS_IJSS2017_Implicit_BE integrator
{
  
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;

  std::ostringstream outs;
  std::string localstr="";
  
  // Model integrator test
  double Dt = 0.01;
  double InitialTime=0, FinalTime = 10000 * Dt;

  // 0
  // Set the time-frame
  this->SetTimeFrame( InitialTime, FinalTime, Dt, 0 );

  // 1
  // Test for S0
  double CurrentTime = 0;
  
  double p0 = this->Parameters->K_p0 ;
  double HardLawJp0Coeff = pow(exp( this->HardeningLaw(p0) ), 1.0/(double) dim);
  
  FTensors Id2 = ttl::identity( i,j ); // ttl::Delta<2,dim,double>(1);
  FTensors pF0 = Id2( i,j ) * HardLawJp0Coeff;
  FTensors eF0 = F0;
  FTensors S0 = this->SecondPKTensor( eF0, p0 );
  
  std::cout << "Initialization: \n";
  this->Initialize( p0, F0, pF0, S0 );
  this->PrintStep( this->timestep, CurrentTime, true );
  
  
  
  // 2
  // Test for S derivative wrt M at S0
  // DSDM( const FTensors& F, const FTensors& eF, const FTensors& pF, const double pc )
  // F1 is the (arbitrary) tensor Fnp1
  
  {
    std::cout << "\n" << "S derivative wrt M at S0 for Fn+1 = { 0.9325570876823829, 0.09245584650000224, 0.09585799776330175, 0.09101817722825798, 0.9681530821466998, 0.09616855183330242, 0.09477765111130806, 0.09159735741693511, 0.9353780514904113 }: \n";
    
    ttl::Tensor<4, dim, double> dsdm = DSDM( F1, eF0, pF0, p0 );
    
    // scrivere dsdm che e` un tensore del quarto ordine
    for ( unsigned i=0; i<dim; i++ )
      for ( unsigned j=0; j<dim; j++ )
      {
        localstr = "";
        PrintInVectorForm(dsdm, i, j, localstr);
        std::cout << " dsdm test for dimension " << dim << ". i =" << i << ", j =" << j << ", dsdm(i,j) = " << localstr;
      }
  }

  
  // 3
  // Test for Gamma derivative wrt M at S0
  // DGammaDM( const ttl::Tensor<4, dim, double>& DSDM, const FTensors& M, const FTensors& S)
  // F1 is the (arbitrary) tensor Fnp1

  {
    FTensors InvF1 = ttl::inverse( F1 );
    FTensors M0;
    
    {
      static constexpr ttl::Index<'i'> i;
      static constexpr ttl::Index<'j'> j;
      static constexpr ttl::Index<'k'> k;
      static constexpr ttl::Index<'l'> l;
      M0 ( i,j ) = pF0( i,k ) * InvF1 ( k,l ) * eF0( l,j );
    }
    
    std::cout << "\n" << "Gamma derivative wrt M at S0 for \n" ;
    std::string tempstr = "", Mtempstr = "";
    PrintInVectorForm( F1, tempstr);
    PrintInVectorForm( M0, Mtempstr);
    std:: cout << "Fn+1 = " << tempstr << ", and \n";
    std:: cout << "M0 = " << Mtempstr << " : \n";
    
    ttl::Tensor<4, dim, double> dsdm = DSDM( F1, eF0, pF0, p0 );
    ttl::Tensor<4, dim, double> dgammadm = DGammaDM( dsdm, M0, S0 );
    
    // scrivere dsdm che e` un tensore del quarto ordine
    for ( unsigned i=0; i<dim; i++ )
      for ( unsigned j=0; j<dim; j++ )
      {
        localstr = "";
        PrintInVectorForm(dgammadm, i, j, localstr);
        std::cout << " dgammadm test for dimension " << dim << ". i =" << i << ", j =" << j << ", dsdm(i,j) = " << localstr;
      }
  }
 
  
  // 4
  // Test for Kirchoff stress derivative wrt M at S0
  // DKSDM( const ttl::Tensor<4, dim, double>& dgammadm, const FTensors& F1, const FTensors& pF)
  // F1 is the (arbitrary) tensor Fnp1
  
  {
    FTensors InvF1 = ttl::inverse( F1 );
    FTensors M0;
    
    {
      static constexpr ttl::Index<'i'> i;
      static constexpr ttl::Index<'j'> j;
      static constexpr ttl::Index<'k'> k;
      static constexpr ttl::Index<'l'> l;
      M0 ( i,j ) = pF0( i,k ) * InvF1 ( k,l ) * eF0( l,j );
    }
    
    std::cout << "\n" << "KS derivative wrt M at S0 for \n" ;
    std::string tempstr = "", Mtempstr = "";
    PrintInVectorForm( F1, tempstr);
    PrintInVectorForm( M0, Mtempstr);
    std:: cout << "Fn+1 = " << tempstr << ", and \n";
    std:: cout << "M0 = " << Mtempstr << " : \n";
    
    ttl::Tensor<4, dim, double> dsdm = DSDM( F1, eF0, pF0, p0 );
    ttl::Tensor<4, dim, double> dgammadm = DGammaDM( dsdm, M0, S0 );
    ttl::Tensor<4, dim, double> dtaudm = DKSDM( dgammadm, F1, pF0 );
    
    
    // scrivere dsdm che e` un tensore del quarto ordine
    for ( unsigned i=0; i<dim; i++ )
      for ( unsigned j=0; j<dim; j++ )
      {
        localstr = "";
        PrintInVectorForm(dtaudm, i, j, localstr);
        std::cout << " dtaudm test for dimension " << dim << ". i =" << i << ", j =" << j << ", dsdm(i,j) = " << localstr;
      }
  }
 
  
  // 5
  // Test for DnormdevtauDtau
  // ttl::Tensor<2, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DnormdevtauDtau( const FTensors& A, const bool ZeroWarning )
  // KS is the Kirchoff stress tensor at S0
  
  {
    
    static constexpr ttl::Index<'a'> a;
    static constexpr ttl::Index<'b'> b;
    static constexpr ttl::Index<'x'> x;
    static constexpr ttl::Index<'z'> z;

    FTensors InvF1 = ttl::inverse( F1 );
    FTensors M0;
    
    {
      static constexpr ttl::Index<'i'> i;
      static constexpr ttl::Index<'j'> j;
      static constexpr ttl::Index<'k'> k;
      static constexpr ttl::Index<'l'> l;
      M0 ( i,j ) = pF0( i,k ) * InvF1 ( k,l ) * eF0( l,j );
    }

    // Kirchoff stress at iteration r
    FTensors KSr;
    KSr( a,b ) = det( pF0 ) / det( M0 ) * eF0( a,x ) * S0( x,z ) * eF0( b,z );
    FTensors DtbarDtau = DnormdevtauDtau( KSr );
  
    std::cout << "\n" << "DnormdevtauDtau at S0 for " ;

    localstr = "";
    PrintInVectorForm(KSr, localstr);
    std::cout << " Kirchoff Stress = "   << localstr;

    localstr = "";
    PrintInVectorForm(DtbarDtau, localstr );
    std::cout << " DnormdevtauDtau = "   << localstr;
  
  }
  
  // 6
  // Test for DgammadotdDM at S0
  // ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DgammadotdDM(  const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& Mr, const double pcr )
  // F1 is the (arbitrary) tensor Fnp1
  
  {
    
    FTensors InvF1 = ttl::inverse( F1 );
    FTensors M0;
    
    {
      static constexpr ttl::Index<'i'> i;
      static constexpr ttl::Index<'j'> j;
      static constexpr ttl::Index<'k'> k;
      static constexpr ttl::Index<'l'> l;
      M0 ( i,j ) = pF0( i,k ) * InvF1 ( k,l ) * eF0( l,j );
    }
    
    std::cout << "\n" << "DgammadotdDM at S0 for " ;
    std::string tempstr = "", Mtempstr = "";
    PrintInVectorForm( F1, tempstr);
    PrintInVectorForm( M0, Mtempstr);
    std:: cout << "Fn+1 = " << tempstr << ", and M0 = " << Mtempstr;

    ttl::Tensor<2, dim, double> dgmdotdm = DgammadotdDM( pF0, F1, eF0, S0, M0, p0 );
    
    localstr = "";
    PrintInVectorForm(dgmdotdm, localstr);
    std::cout << " DgammadotdDM test for dimension " << dim  << localstr;
  }

  // 7
  // Test for DPsidDM at S0
  // ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DPsidDM( const FTensors& devS, const FTensors& F, const FTensors& eF, const FTensors& pF, const double pc )
  // F1 is the (arbitrary) tensor Fnp1
  
  {
    
    // deviatoric part of the Second-Piola Kirchoff stress
    FTensors devS = dev( S0 );
    
    std::cout << "\n" << "DPsidDM at S0 for " ;
    std::string tempstr = "", Mtempstr = "";
    PrintInVectorForm( F1, tempstr);
    std:: cout << "Fn+1 = " << tempstr;
    
    ttl::Tensor<4, dim, double> dpsidm = DPsidDM( devS, F1, eF0, pF0, p0 );
    
    // scrivere dsdm che e` un tensore del quarto ordine
    for ( unsigned i=0; i<dim; i++ )
      for ( unsigned j=0; j<dim; j++ )
      {
        localstr = "";
        PrintInVectorForm(dpsidm, i, j, localstr);
        std::cout << " DPsidDM test for dimension " << dim << ". i =" << i << ", j =" << j << ", dsdm(i,j) = " << localstr;
      }
  }
  
  
  
  // 8
  
  // Test for dRMdM at S0
  // ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::dRMdM( const double Deltat, const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& Mr, const double pcr )
  // F1 is the (arbitrary) tensor Fnp1
  
  {
    
    FTensors InvF1 = ttl::inverse( F1 );
    FTensors M0;
    
    {
      static constexpr ttl::Index<'i'> i;
      static constexpr ttl::Index<'j'> j;
      static constexpr ttl::Index<'k'> k;
      static constexpr ttl::Index<'l'> l;
      M0 ( i,j ) = pF0( i,k ) * InvF1 ( k,l ) * eF0( l,j );
    }
    
    std::cout << "\n" << "dRMdM at S0 for " ;
    std::string tempstr = "", Mtempstr = "";
    PrintInVectorForm( F1, tempstr);
    PrintInVectorForm( M0, Mtempstr);
    std:: cout << "Fn+1 = " << tempstr << ", and M0 = " << Mtempstr;
    
    ttl::Tensor<4, dim, double> drdm = DRMDM( this->Deltat, pF0, F1, eF0, S0, M0, p0 );
    
    // scrivere dsdm che e` un tensore del quarto ordine
    for ( unsigned i=0; i<dim; i++ )
      for ( unsigned j=0; j<dim; j++ )
      {
        localstr = "";
        PrintInVectorForm(drdm, i, j, localstr);
        std::cout << " DPsidDM test for dimension " << dim << ". i =" << i << ", j =" << j << ", drdm(i,j) = " << localstr;
      }
  }

  
  // 9
  
  // Test for DRpcDpc at S0
  // double KMS_IJSS2017_Implicit_BE<dim>::DRpcDpc( const double Deltat, const FTensors& pFn, const FTensors& Fnp1, const FTensors& eFr, const FTensors& Sr, const FTensors& Mr, const double pcr )
  // F1 is the (arbitrary) tensor Fnp1
  
  {
    
    FTensors InvF1 = ttl::inverse( F1 );
    FTensors M0;
    
    {
      static constexpr ttl::Index<'i'> i;
      static constexpr ttl::Index<'j'> j;
      static constexpr ttl::Index<'k'> k;
      static constexpr ttl::Index<'l'> l;
      M0 ( i,j ) = pF0( i,k ) * InvF1 ( k,l ) * eF0( l,j );
    }
    
    std::cout << "\n" << "dRMdM at S0 for " ;
    std::string tempstr = "", Mtempstr = "";
    PrintInVectorForm( F1, tempstr);
    PrintInVectorForm( M0, Mtempstr);
    std:: cout << "Fn+1 = " << tempstr << ", and M0 = " << Mtempstr;
    
    double drdpc = DRpcDpc( this->Deltat, pF0, F1, eF0, S0, M0, p0 );
    
    // scrivere drdpc 
    std::cout << " drdpc test for dimension " << dim << ". drdpc =" << drdpc << "\n";
  }
  


}



template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017_Implicit_BE<dim>::UpdateHattedZwx( FTensors M )
//! Updates operator Zwx
//! See page 17b of the "Consistent Tangent Operator" notes
//! See page 1 of the "Debugging and code optimization" notes
{
  
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'x'> x;
    
  // Hardening function at this->pcnp1
  double Hpc = this->HardeningLaw( this->pcnp1 );
  
  // Inverse of the determinant
  double detMm1 = 1.0 / det(M);
  
  // In matrix calculus, Jacobi's formula expresses the derivative of the determinant of a matrix A in terms of the adjugate of A.
  // d Det[A] / dA_ij = adjT(A)_ij
  // the tensor adjTM is the transpose of the adjugate of M
  ttl::Tensor<2, dim, double> adjTM =  adjugateTranspose( M );
  
  double ZwxCoeff = det( this->pFn ) * this->pcnp1 * this->pcnp1 * detMm1 * detMm1 /  exp( Hpc ) ;
  
  FTensors Zwx =  ZwxCoeff * adjTM( w,x );
  return Zwx;
  
  
}

template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017_Implicit_BE<dim>::UpdateZwx( FTensors M )
//! Updates operator Zwx
//! See page 17b of the Consistent Tangent Operator notes
{
  
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'x'> x;
  
  // Material parameters
  double a1=this->Parameters->hr_a1, a2=this->Parameters->hr_a2, l1=this->Parameters->hr_Lambda1, l2=this->Parameters->hr_Lambda2;
  
  // Hardening function at this->pcnp1
  double Hpc = this->HardeningLaw( this->pcnp1 );
  
  // derivative dHpc/dpc
  double dHpcdpc = a1 * l1 * exp( - l1 / this->pcnp1 ) + a2 * l2 * exp( - l2 / this->pcnp1 ) ;
  
  // Inverse of the determinant
  double detMm1 = 1.0 / det(M);
  
  // In matrix calculus, Jacobi's formula expresses the derivative of the determinant of a matrix A in terms of the adjugate of A.
  // d Det[A] / dA_ij = adjT(A)_ij
  // the tensor adjTM is the transpose of the adjugate of M
  ttl::Tensor<2, dim, double> adjTM =  adjugateTranspose( M );
  
  double ZwxCoeff = det( this->pFn ) * this->pcnp1 * this->pcnp1 * detMm1 * detMm1 / (  dHpcdpc *  exp( Hpc ) ) ;

  FTensors Zwx =  ZwxCoeff * adjTM( w,x );
  return Zwx;

  
}




template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DSDpc()
//! Derivative of S wrt pc
//! see page 23 of the Consistent tangent operator notes
{
  
  using namespace ttlindexes;
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'y'> y;

  FTensors DSisoDpc, DSvolDpc;
  
  // S isochoric
  
  double Je = det( this->eFnp1 ), onethird = 1.0/3.0, powJe =  pow(Je, - 2.0 * onethird ) ;
  
  // Right Cauchy-Green tensor ...
  FTensors Ce, InvCe;
  Ce( z,y ) = this->eFnp1( k,z ) * this->eFnp1( k,y );
  InvCe = ttl::inverse( Ce );
  
  DSisoDpc( z,y ) = this->DmuDpc( this->pcnp1 ) * powJe * ( ttl::identity<dim>( z,y ) - onethird * Trace( Ce ) * InvCe( z,y )  );

  // S volumetric

  double coh = this->c( this->pcnp1 );
  double alpha = pow( 1.0 -  coh / this->Parameters->c_inf, this->Parameters->pl_n ) / this->Parameters->K_kappa  ;

  DSvolDpc( z,y ) = (
                     0.5 * this->DbulkDpc( this->pcnp1 ) * ( Je * ( log(Je) + 1 ) - 1 ) +
                     this->DcDpc( this->pcnp1 ) * ( 1 - pow( Je, -alpha ) ) + ( coh + this->Parameters->K_p0 ) * alpha * pow( Je, -alpha-1 ) * this->DalphaDpc( this->pcnp1 )
                    ) * InvCe( z,y ) ;
  
  return DSisoDpc( z,y ) + DSvolDpc( z,y );

}



template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::dSdF( const FTensors& F, const FTensors& eF, const FTensors& pF, const double pc )
//! Derivative of S wrt F
//! see page 25 of the Consistent tangent operator notes
//! This derivative defines the 4th order tensor V(zykl)
//!
//! Note that this is the partial derivative of S wrt F, the total derivative
//! accounts also for DSDM DMDF.
//!
{
 
  using namespace ttlindexes;
  
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'e'> e;
  static constexpr ttl::Index<'f'> f;
  static constexpr ttl::Index<'g'> g;
  static constexpr ttl::Index<'h'> h;
  static constexpr ttl::Index<'x'> x;
  
  
  FTensors pFinv;
  pFinv = ttl::inverse( pF );
  
  double Je = det( eF ), onethird = 1.0/3.0, powJe =  pow(Je, - 2.0 * onethird ) ;
  double coh = this->c( pc );
  double alpha = pow( 1.0 -  coh / this->Parameters->c_inf, this->Parameters->pl_n ) / this->Parameters->K_kappa  ;
  
  // delta
  ttl::Tensor<4, dim, double> Id4 = ttl::identity( i,j,k,l ); // ttlIdentity<dim,double>(); // ttl::Delta<4,dim,double>(1); //
  ttl::Tensor<2, dim, double> Id2 = ttl::identity( i,j ); // ttl::Delta<2,dim,double>(1);
  
  // Left Caucgy-Green tensor ...
  FTensors Ce, InvCe;
  Ce(e,f) = eF(x,e) * eF(x,f);
  InvCe = ttl::inverse( Ce );
  
  // ... and its derivative
  ttl::Tensor<4, dim, double> DCeDeF;
  DCeDeF(e,f,g,h) = Id4(x,e,g,h) * eF(x,f) + eF(x,e) * Id4(x,f,g,h);
  
  
  // derivatives of the Second Piola-Kirchoff stress wrt Ce
  ttl::Tensor<4, dim, double> DSisoDCe;
  DSisoDCe(a,b,e,f) =
  this->shearmodulus( pc ) * ( - onethird * powJe *  InvCe(e,f) ) * ( Id2(a,b) - Trace( Ce ) * onethird * InvCe(a,b) ) +
  this->shearmodulus( pc ) * powJe * onethird * ( -  Id2(e,f) * InvCe(a,b) + Trace( Ce ) * InvCe(a,e) * InvCe(f,b) );
  
  ttl::Tensor<4, dim, double> DSvolDCe;
  DSvolDCe(a,b,e,f) =
  ( this->bulkmodulus( pc ) + alpha * pow( Je, -alpha-1 ) * ( coh + this->Parameters->K_p0 ) + 0.5 * this->bulkmodulus( pc ) * log( Je ) ) * 0.5 * Je * InvCe(e,f) * InvCe(a,b)
  -( coh + 0.5 *  this->bulkmodulus( pc ) * ( Je -1 ) - pow( Je, -alpha ) * ( coh + this->Parameters->K_p0 ) + 0.5 * this->bulkmodulus( pc ) * Je * log( Je )  ) * InvCe(a,e) * InvCe(f,b);

  ttl::Tensor<4, dim, double> DSDCe;
  DSDCe(a,b,e,f) = DSisoDCe(a,b,e,f) + DSvolDCe(a,b,e,f);
  
  // derivative of the eF tensor wrt F
  ttl::Tensor<4, dim, double> DeFDF;
  DeFDF( g,h,k,l ) = Id4( g,x,k,l ) * pFinv( x,h );
  
  // DSDF
  ttl::Tensor<4, dim, double> buf;
  buf( a,b,k,l ) =  DSDCe( a,b,e,f ) * DCeDeF( e,f,g,h ) * DeFDF( g,h,k,l );
  
  return buf;


}


template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::Minor_cdkl(
                                                                      const FTensors& F,
                                                                      const FTensors& eF,
                                                                      const FTensors& pF,
                                                                      const FTensors& S,
                                                                      const FTensors& M,
                                                                      ttl::Tensor<4, dim, double> Vzykl
                                                                    )
//! Second order tensor Minor_cdkl()
//! see page 17 of the Consistent tangent operator notes
{
  
  using namespace ttlindexes;
  static constexpr ttl::Index<'c'> c;
  static constexpr ttl::Index<'d'> d;
  static constexpr ttl::Index<'e'> e;
  static constexpr ttl::Index<'f'> f;
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'y'> y;
  
  // definition of tensors Aa and Gg, based on the definition of the Kirchoff
  // stress as at page 1 of the Constitutive staggered integrator notes
  // the quantity det( this->pFn ) has been taken out of the definition of Aa
  // Note that Bb has not been defined since it is the transpose of Aa
  
  FTensors Gg;
  double detMm1 = 1.0 / det(M);
  Gg( z,w ) =  M( z,i ) * S( i,j ) * M( w,j ) * detMm1 ;
  
  FTensors Aa, InvFp = ttl::inverse( pF );
  Aa( i,z ) =  F( i,x ) * InvFp( x,z );

  // auxiliary tensors
  FTensors InvpFn = ttl::inverse( this->pFn );

  // delta
  ttl::Tensor<4, dim, double> Id4 = ttl::identity( i,j,k,l );
  
  ttl::Tensor<4, dim, double> Minor = {};
  
  Minor( c,d,k,l ) = det( this->pFn ) *
                    (
                     Id4( c,z,k,l ) * InvpFn( z,e ) * Gg( e,f ) * Aa( d,f ) +
                     Aa( c,e ) * Gg( e,f ) * InvpFn( z,f ) * Id4( d,z,k,l ) +
                     Aa( c,e ) * detMm1 * M( e,z ) * Vzykl( z,y,k,l ) * M( f,y ) * Aa( d,f )
                    );
  
  return Minor;

  
}



template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::HattedMajor_cdwx(
                                                                      const FTensors& F,
                                                                      const FTensors& eF,
                                                                      const FTensors& pF,
                                                                      const FTensors& S,
                                                                      const FTensors& M,
                                                                      ttl::Tensor<4, dim, double> HattedLzywx,
                                                                      const double dHpdpc
                                                                      )
//! Second order tensor HattedMajor_cdwx
//! See page 17b of the Consistent tangent operator notes
//! See page 2 of the "Debugging aand code optimization" notes
{
  
  // we take advantage of the functions implemented for the integrator.
  // note that dsdm has been replaced by Lzywx, which incorporates DSDpc times Zwx
  
  ttl::Tensor<4, dim, double> dgammadm = HattedDGammaDM( HattedLzywx, dHpdpc, M, S );
  ttl::Tensor<4, dim, double> dtaudm = DKSDM( dgammadm, F, pF );
  
  return dtaudm;
  
}




template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::Major_cdwx(
                                                                      const FTensors& F,
                                                                      const FTensors& eF,
                                                                      const FTensors& pF,
                                                                      const FTensors& S,
                                                                      const FTensors& M,
                                                                      ttl::Tensor<4, dim, double> Lzywx
                                                                      )
//! Second order tensor Major_cdwx()
//! see page 17b of the Consistent tangent operator notes
{
  
  using namespace ttlindexes;
  static constexpr ttl::Index<'c'> c;
  static constexpr ttl::Index<'d'> d;
  static constexpr ttl::Index<'e'> e;
  static constexpr ttl::Index<'f'> f;
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'y'> y;
  
  // we take advantage of the functions implemented for the integrator.
  // note that dsdm has been replaced by Lzywx, which incorporates DSDpc times Zwx
  
  ttl::Tensor<4, dim, double> dgammadm = DGammaDM( Lzywx, M, S );
  ttl::Tensor<4, dim, double> dtaudm = DKSDM( dgammadm, F, pF );

  return dtaudm;

}


template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::Phi_abzy( const FTensors& S, bool ZeroWarning )
//! Fourth order tensor Phi_abzy()
//! see page 18 of the Consistent tangent operator notes
{
  static constexpr ttl::Index<'a'> a;
  static constexpr ttl::Index<'b'> b;
  static constexpr ttl::Index<'c'> c;
  static constexpr ttl::Index<'d'> d;
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'y'> y;
  
  
  // delta
  FTensors Id2 = ttl::identity<dim>( a,b );

  FTensors devS;
  devS( z,y ) = dev( S )( z,y );
  
  double normdevS = FrobeniusNorm( devS );
  
  if ( std::fpclassify(normdevS) == FP_ZERO )
  {
    if ( ZeroWarning )
    {
      std::cerr << " WARNING: KMS_IJSS2017_Implicit_BE<dim>::Phi_abzy( const FTensors& S) ";
      std::cerr << " Derivative of the deviatoric part not defined since || dev( S ) || = " << normdevS << ". \n";
      std::cerr << " Code did not abort, returned ZERO but outcomes might be wrong. \n";
    }
      
    return ttl::zero<dim>( a,b,c,d );
  }
      
  double den = 1.0 / normdevS;
      
  ttl::Tensor<4, dim, double> Phi;
  
  Phi( a,b,z,y ) = den * ( Id2( a,c ) * Id2( b,d ) - den * den * devS( a,b ) * devS( c,d ) ) * ( Id2( c,z ) * Id2( d,y ) - Id2( c,d ) * Id2( z,y ) / 3.0 );

  return Phi;
      
}






template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE<dim>::DSDF(
                                                                ttl::Tensor<4, dim, double> dmdf,
                                                                const double dt,
                                                                const bool Verbose
                                                                )
//! This function evaluates the fourth order tensor DS/DF (total derivative) with dmdf passed as an argument
//! This method expects that the Step has already been updated elsewhere, since it uses quantities at step n+1
//! that MUST have been updated before evaluating dmdf
//!
//! The optimal method would evaluate DSDF and DMDF in the same function, since the two functions share
//! several common tensors. This has been implemented in the method DMDFandDSDF
//!
//! 1 - Calculate intermediate operators
//! 2 - Return
{
  
  using namespace ttlindexes;
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'y'> y;
  
  
  // 1 - Calculate intermediate operators
  
  // tensor M = pFn * InvpFnp1
  FTensors InvpFnp1 = ttl::inverse( this->pFnp1 );
  FTensors M = this->pFn( i,k ) * InvpFnp1( k,j );
  
  // Zwx has been defined at page 17b of the consistent tangent operator notes
  FTensors Zwx = UpdateZwx( M );
  
  // Lambdazywx has been defined at page 16 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> Lambdazywx;
  ttl::Tensor<4, dim, double> dsdm = DSDM( this->Fnp1, this->eFnp1, this->pFnp1, this->pcnp1 );
  Lambdazywx( z,y,w,x ) = dsdm( z,y,w,x ) + DSDpc() ( z,y ) * Zwx( w,x );
  
  // Vzykl has been defined at page 16 of the consistent tangent operator notes
  // Note that dSdF is the partial derivative and differs from DSDF
  ttl::Tensor<4, dim, double> Vzykl = dSdF( this->Fnp1, this->eFnp1, this->pFnp1, this->pcnp1 );
  
  
  // 2 - Return
  
  ttl::Tensor<4, dim, double> dsdf;
  dsdf( z,y,k,l ) = Lambdazywx( z,y,w,x ) * dmdf( w,x,k,l ) + Vzykl( z,y,k,l );
  
  return dsdf;
  
}

















#endif /* KMS_IJSS2017_implicit_templates_hpp */
