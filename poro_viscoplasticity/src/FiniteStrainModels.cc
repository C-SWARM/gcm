//
//  FiniteStrainModels.cpp
//  ttl-learning
//
//  Created by Alberto Salvadori on 2/2/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#include "FiniteStrainModels.h"



// class FiniteStrainModels_IO
// ***************************


// Constructors
// ------------



FiniteStrainModels_IO::FiniteStrainModels_IO(const std::string& logstr,
                                 const std::string& Fstr,
                                 const std::string& SPKstr,
                                 const std::string& KSstr,
                                 const std::string& sigmastr,
                                 unsigned nts) :
MechanicalModels_IO ( logstr, sigmastr, nts ),
Ffile( Fstr, std::ios::out),
SPKfile( SPKstr, std::ios::out),
KSfile( KSstr, std::ios::out)
{}



FiniteStrainModels_IO::~FiniteStrainModels_IO()
{
  Ffile.close();
  SPKfile.close();
  KSfile.close();
}
