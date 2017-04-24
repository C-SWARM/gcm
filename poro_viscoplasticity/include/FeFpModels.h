//
//  FeFpModels.hpp
//  ttl-learning
//
//  Created by Alberto Salvadori on 2/2/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#ifndef FeFpModels_hpp
#define FeFpModels_hpp

#include <stdio.h>
#include "FiniteStrainModels.h"


class FeFpModels_IO : public FiniteStrainModels_IO
/*!
 \brief This class contains and handles the input and output for the KMS-IJSS2017 model
 */
{
protected:
  
  std::ofstream pFfile;     //!< output of Fp for FeFpModels
  
public:
  // Constructors & destructors
  FeFpModels_IO(){}
  FeFpModels_IO( const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, unsigned = 1 );
  virtual ~FeFpModels_IO();
  
  // Accessing IO
  std::ofstream& pF(){ return pFfile; }
  
};



#endif /* FeFpModels_hpp */
