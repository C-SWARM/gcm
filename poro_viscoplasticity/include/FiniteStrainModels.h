//
//  FiniteStrainModels.h
//  ttl-learning
//
//  Created by Alberto Salvadori on 2/2/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#ifndef FiniteStrainModels_h
#define FiniteStrainModels_h

#include <stdio.h>
#include "MechanicalModels.h"



class FiniteStrainModels_IO : public MechanicalModels_IO
/*!
 \brief This class contains and handles the input and output for the KMS-IJSS2017 model
 */
{
protected:
  std::ofstream Ffile;      //!< output of F for KMS_IJSS2017
  std::ofstream SPKfile;    //!< output of the Second Piola-Kirchoff stress for KMS_IJSS2017
  std::ofstream KSfile;     //!< output of the Kirchoff stress for KMS_IJSS2017
  
  
public:
  // Constructors & destructors
  FiniteStrainModels_IO(){}
  FiniteStrainModels_IO( const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, unsigned = 1 );
  virtual ~FiniteStrainModels_IO();
  
  // Accessing IO
  std::ofstream& F(){ return Ffile; }
  std::ofstream& SPK(){ return SPKfile; }
  std::ofstream& KS(){ return KSfile; }
  
};





#endif /* FiniteStrainModels_hpp */
