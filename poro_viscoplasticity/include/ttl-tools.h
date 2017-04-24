//
//  KMS-IJSS2017-tools.h
//  ttl-learning
//
//  Created by Alberto Salvadori on 1/16/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

/*!
 
 \file ttl-tools.h
 \brief Header file for the ttl-tools model
 
 This file contains some ttl-methods that are necessary for the integration of the mechanical models.
 they are:
 - Trace operator
 - Dev operator
 - Frobenius norm
 - Adjugate
 
 */


#ifndef ttl_tools_h
#define ttl_tools_h

#include <stdio.h>
#include <iostream>
#include <iomanip>

#include <math.h>

namespace ttlindexes
{
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;
}


// outputs 2nd order tensors in vector form
template <int R, int D, class S = double>
void
PrintInVectorForm(const ttl::Tensor<R,D,S>& A, std::string& out)
{
  
  std::ostringstream outs;
  for (int r=0;r<pow(D,R)-1;r++)
    outs << std::setw(20) << std::setprecision(15) << A.get(r) << ", ";
  outs << std::setw(20) << std::setprecision(15) << A.get(pow(D,R)-1) << "\n";
  out +=  outs.str();
  
}


// outputs 2nd order tensors in vector form
template <int R, int D, class S = double>
void
PrintInMathematicaForm(const ttl::Tensor<R,D,S>& A, std::string& out)
{
  
  std::ostringstream outs;
  
  outs << "{";
  int r=0;
  for (int i =0; i<pow(D,R/2); i++)
  {
    outs << "{";
    for (int j =0; j<pow(D,R/2)-1; j++)
      outs << std::setw(20) << std::setprecision(15) << A.get( r++ ) << ", ";
    
    outs << std::setw(20) << std::setprecision(15) << A.get( r++ ) << " ";
    
    if ( i == pow(D,R/2) )
      outs <<  "}";
    else
      outs << "},";
  }
  outs <<  "}";

  out +=  outs.str();
  
}


template <int dim, class S = double>
void
PrintInVectorForm(const ttl::Tensor<4, dim, S>& A, unsigned i, unsigned j, std::string& out)
//! This method prints the NxN subtensor i,j of the fourth order tensor A
{
  
  std::ostringstream outs;
  unsigned init = ( i * dim  + j ) * dim * dim ;
  for (int r=0; r< dim * dim - 1 ; r++)
      outs << std::setw(20) << std::setprecision(15) << A.get( init + r ) << ", ";
  outs << std::setw(20) << std::setprecision(15) << A.get( init + pow(dim, 2)-1) << "\n";
  out +=  outs.str();
  
}


template <int dim>
double Trace( const ttl::Tensor<2, dim, double> T )
{
  static constexpr ttl::Index<'i'> i;
  return T(i,i);
}



template <int dim>
ttl::Tensor<2, dim, double> dev( const ttl::Tensor<2, dim, double> T )
{
  
  using namespace ttlindexes;
  
  ttl::Tensor<2, dim, double> Id2 = ttl::identity(i,j);

  return T(i,j) - Id2(i,j) * ( Trace(T) / 3.0 );
}



template <int dim>
ttl::Tensor<2, 3, double> dev3D( const ttl::Tensor<2, dim, double> T )
  //! The deviator operator is quite tricky, because for plane strain or plane stress
  //! problems, the deviator does not lay in the same two-dimensional space of the
  //! tensor that generates it. To this purpose, I made the function dev3D, which
  //! returns a three-dimensional tensor independent upon the generator.
  //!
  //! As a typical example, this function must be used in the three-field formulation,
  //! for the deviatoric subspace
{
  
  
  using namespace ttlindexes;
  ttl::Tensor<2, 3, double> Id2 = ttl::identity<3>(i,j);
  
  if ( dim == 3 )
  {
    ttl::Tensor<2, 3, double> T3D = { T(0,0), T(0,1), T(0,2), T(1,0), T(1,1), T(1,2), T(2,0), T(2,1), T(2,2)};
    
    return T3D(i,j) - Id2(i,j) * ( Trace(T3D) / 3.0 );
  }
  
  else if ( dim == 2 )
  {
    ttl::Tensor<2, 3, double> T3D = { T(0,0), T(0,1), 0, T(1,0), T(1,1), 0, 0, 0, 0 };
    return T3D(i,j) - Id2(i,j) * ( Trace(T3D) / 3.0 );
  }
  
  else if ( dim == 1 )
  {
    ttl::Tensor<2, 3, double> T3D = { T(0,0),0, 0, 0, 0, 0, 0, 0, 0 };
    return T3D(i,j) - Id2(i,j) * ( Trace(T3D) / 3.0 );
  };
  
  std::cout << " WARNING: ttl::Tensor<2, 3, double> dev3D( const ttl::Tensor<2, dim, double> T ) - The dimension is not defined.";
  std::cout << " Code did not abort but outcomes might be wrong. \n";
  
  return Id2;

}





template <int dim>
ttl::Tensor<2, 3, double> dev3D( const ttl::Tensor<2, 3, double> T )
//! The deviator operator is quite tricky, because for plane strain or plane stress
//! problems, the deviator does not lay in the same two-dimensional space of the
//! tensor that generates it. To this purpose, I made the function dev3D, which
//! returns a three-dimensional tensor independent upon the generator.
//!
//! As a typical example, this function must be used in the three-field formulation,
//! for the deviatoric subspace
{
  using namespace ttlindexes;
  ttl::Tensor<2, 3, double> Id2 = ttl::identity(i,j);

  return T(i,j) - Id2(i,j) * ( Trace(T) / 3.0 );
}



template <int dim>
ttl::Tensor<2, 3, double> dev3D( const ttl::Tensor<2, 2, double> T )
//! The deviator operator is quite tricky, because for plane strain or plane stress
//! problems, the deviator does not lay in the same two-dimensional space of the
//! tensor that generates it. To this purpose, I made the function dev3D, which
//! returns a three-dimensional tensor independent upon the generator.
//!
//! As a typical example, this function must be used in the three-field formulation,
//! for the deviatoric subspace
{
  using namespace ttlindexes;
  ttl::Tensor<2, 3, double> Id2 = ttl::identity(i,j);
  
  ttl::Tensor<2, 3, double> T3D = { T(0,0), T(0,1), 0, T(1,0), T(1,1), 0, 0, 0, 0 };
  return T3D(i,j) - Id2(i,j) * ( Trace(T) / 3.0 );
}



template <int dim>
ttl::Tensor<2, 3, double> dev3D( const ttl::Tensor<2, 1, double> T )
//! The deviator operator is quite tricky, because for plane strain or plane stress
//! problems, the deviator does not lay in the same two-dimensional space of the
//! tensor that generates it. To this purpose, I made the function dev3D, which
//! returns a three-dimensional tensor independent upon the generator.
//!
//! As a typical example, this function must be used in the three-field formulation,
//! for the deviatoric subspace
{
  using namespace ttlindexes;
  ttl::Tensor<2, 3, double> Id2 = ttl::identity(i,j);
  
  ttl::Tensor<2, 3, double> T3D = { T(0,0), 0, 0, 0, 0, 0, 0, 0, 0 };
  return T3D(i,j) - Id2(i,j) * ( Trace(T) / 3.0 );
}



template <int dim>
double FrobeniusNorm(const ttl::Tensor<2, dim, double> T)
{
  using namespace ttlindexes;
  
  return sqrt( T(i,j)*T(i,j) );
}


template <int dim>
double FrobeniusNorm(const ttl::Tensor<4, dim, double> T)
{
  using namespace ttlindexes;
  
  return sqrt( T( i,j,k,l )*T( i,j,k,l ) );
}


template <int dim>
ttl::Tensor<2, dim, double> adjugateTranspose(const ttl::Tensor<2, dim, double> M)
//! In matrix calculus, Jacobi's formula expresses the derivative of the determinant of a matrix A in terms of the adjugate of A.
//! d Det[A] / dA_ij = adjT(A)_ij
//! the tensor adjTM is the transpose of the adjugate of M

{
  if ( dim==2 )
    return { M.get(3), -M.get(2), -M.get(1), M.get(0) };
  
  if ( dim==3 )
    return { -M.get(5) * M.get(7) + M.get(4) * M.get(8),  M.get(5) * M.get(6) - M.get(3) * M.get(8), -M.get(4) * M.get(6) + M.get(3) * M.get(7),
              M.get(2) * M.get(7) - M.get(1) * M.get(8), -M.get(2) * M.get(6) + M.get(0) * M.get(8),  M.get(1) * M.get(6) - M.get(0) * M.get(7),
             -M.get(2) * M.get(4) + M.get(1) * M.get(5),  M.get(2) * M.get(3) - M.get(0) * M.get(5), -M.get(1) * M.get(3) + M.get(0) * M.get(4) };
  
  std::cout << " WARNING: ttl::Tensor<2, dim, double> adjugateTranspose(const ttl::Tensor<2, dim, double> M) - The dimension is not defined.";
  std::cout << " Code did not abort but outcomes might be wrong. \n";

}







#endif /* ttl_tools_h */
