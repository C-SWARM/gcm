//
//  FeFpModels.cpp
//  ttl-learning
//
//  Created by Alberto Salvadori on 2/2/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#include "FeFpModels.h"




// class FiniteStrainModels_IO
// ***************************


// Constructors
// ------------



FeFpModels_IO::FeFpModels_IO(const std::string& logstr,
                             const std::string& Fstr,
                             const std::string& pFstr,
                             const std::string& SPKstr,
                             const std::string& KSstr,
                             const std::string& sigmastr,
                             unsigned nts) :
FiniteStrainModels_IO ( logstr, Fstr, SPKstr, KSstr, sigmastr, nts ),
pFfile( pFstr, std::ios::out)
{}



FeFpModels_IO::~FeFpModels_IO()
{
  pFfile.close();
}


