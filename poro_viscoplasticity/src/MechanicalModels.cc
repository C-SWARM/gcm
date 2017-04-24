//
//  MechanicalModels.cpp
//  ttl-learning
//
//  Created by Alberto Salvadori on 2/2/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//


/*! \file MechanicalModels.cpp
 \brief
 
 This file contains abstract classes and methods for a generic material model
 
 */




// include

#include <stdio.h>
#include <math.h>

#include "MechanicalModels.h"




// class MechanicalModels_IO
// *************************


// Constructors
// ------------



MechanicalModels_IO::MechanicalModels_IO(const std::string& logstr,
                                 const std::string& sigmastr,
                                 unsigned nts) :
logfile( logstr, std::ios::out),
sigmafile( sigmastr, std::ios::out)
{
  PrintEveryNSteps = nts;
}



MechanicalModels_IO::~MechanicalModels_IO()
{
  logfile.close();
  sigmafile.close();
}


// Methods
// ------------


bool MechanicalModels_IO::StepToPrint( size_t step )
{
  return ( ( step % PrintEveryNSteps ) == 0 );
}




