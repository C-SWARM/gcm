#ifndef H__H__CONSTITUTIVE_MODEL__H__H
#define H__H__CONSTITUTIVE_MODEL__H__H

#define DEBUG_PRINT_STAT 0

#define DIM_3       3
#define DIM_3x3     9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81

#include <ttl/ttl.h>
#include <math.h>
  
namespace {
  template<int R, int D = DIM_3, class S = double>
  using Tensor  = ttl::Tensor<R, D, S>;
  
  template<int R, int D = DIM_3, class S = double*>
  using TensorA = ttl::Tensor<R, D, S>;

  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'m'> m;
  static constexpr ttl::Index<'r'> r;
  static constexpr ttl::Index<'s'> s;
  static constexpr ttl::Index<'v'> v;
  static constexpr ttl::Index<'w'> w;        
    
  template<class T1, class T2> inline int inv(T1 A, T2 AI)
  {
    int err = 0;

    try
    {
      // attempt to take the inverse
      AI(i,j) = ttl::inverse(A)(i,j);
    }
    catch(const int inverseException)
    {
      // no inverse exists
      err += 1;
    }
    return err;
  }
      
}
    
#endif
