//
//  KMS-IJSS2017.cpp
//  ttl-learning
//
//  Created by alberto salvadori on 12/13/16.
//  Copyright © 2016 alberto salvadori. All rights reserved.
//

/*! \file KMS-IJSS2017.cpp
 \brief A. Krairi, K. Matouš, A. Salvadori, A poro-viscoplastic constitutive model for granular materials at finite strain, submitted to IJSS (2017)
 
 This file contains classes and methods for the poroviscoplastic model desccribed in the paper
 "A. Krairi, K. Matouš, A. Salvadori, A poro-viscoplastic constitutive model for granular materials at finite strain, submitted to IJSS (2017)".
 */


// include
#include <math.h>

#include "KMS-IJSS2017.h"

// class KMS_IJSS2017_IO
// *********************


// Constructors
// ------------



KMS_IJSS2017_IO::KMS_IJSS2017_IO(const std::string& logstr,
                                 const std::string& Fstr,
                                 const std::string& pFstr,
                                 const std::string& SPKstr,
                                 const std::string& KSstr,
                                 const std::string& sigmastr,
                                 const std::string& pcstr,
                                 const std::string& matstr,
                                 unsigned nts) :
FeFpModels_IO( logstr, Fstr, pFstr, SPKstr, KSstr, sigmastr, nts ),
pcfile( pcstr, std::ios::out),
matfile( matstr, std::ios::out)
{}



KMS_IJSS2017_IO::~KMS_IJSS2017_IO()
{
  pcfile.close();
  matfile.close();
}


























