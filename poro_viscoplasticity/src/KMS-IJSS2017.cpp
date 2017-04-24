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



// class KMS_IJSS2017_Parameters
// *****************************


// Constructors
// ------------

KMS_IJSS2017_Parameters::KMS_IJSS2017_Parameters
(double M, double alpha, double m, double gamma_0, double a1, double a2, double L1, double L2,
 double cinf, double Gamma, double B, double pcb, double mu0, double mu1, double p0, double kappa,
 double n, double g0, double pcinf, bool smM)
//! This constructor reads and assigns the material parameters
//! It also assigns the smoothMacaulay bool variable, based on which the functions c and d will
//! use a smoothed Macauley brackets or not.
{
  
  // assigning the material parameters
  yf_M = M;
  yf_alpha = alpha;
  flr_m = m;
  flr_gamma0=gamma_0;
  hr_a1 = a1;
  hr_a2 = a2;
  hr_Lambda1 = L1;
  hr_Lambda2 = L2;
  c_inf = cinf;
  c_Gamma = Gamma;
  d_B = B;
  d_pcb = pcb;
  mu_0 = mu0;
  mu_1 = mu1;
  K_p0 = p0;
  K_kappa = kappa;
  pl_n = n;
  cf_g0 = g0;
  cf_pcinf = pcinf;
  
  // smooth Macauley brackets
  smMbrackets = smM;
}


void KMS_IJSS2017_Parameters::Checks( bool Verbose )
//! This method performs some checks and prints warnings, foreseeing issues in the
//! KMS_IJSS2017 model integration
{
  if ( Verbose && ( hr_a1 * hr_Lambda1 <=0 || hr_a2 * hr_Lambda2 <=0 ) )
  {
    std::cout << " WARNING: KMS_IJSS2017_Explicit<dim>::FindpcFromJp( ) - The sign of the first derivative is not defined.";
    std::cout << " Code did not abort but outcomes might be wrong. \n";
  }
}




// IO - methods
// ------------

void KMS_IJSS2017_Parameters::AsAString( std::string& str)
//! This method prints the model features as a string
{
  
  str += "\n";
  str += "  KMS_IJSS2017 model material parameters:\n";
  str += "   Yield function parameters:\n";
  str += "    M = " + std::to_string(yf_M);
  str += ", alpha = " + std::to_string(yf_alpha);
  str += "\n   Flow rule parameters:\n";
  str += "    m = " + std::to_string(flr_m);
  str += ", gamma0 = " + std::to_string(flr_gamma0);
  str += "\n   Hardening rule parameters:\n";
  str += "    a1 = " + std::to_string(hr_a1);
  str += ", a2 = " + std::to_string(hr_a2);
  str += ", Lambda1 = " + std::to_string(hr_Lambda1);
  str += ", Lambda2 = " + std::to_string(hr_Lambda2);
  str += "\n   Cohesion rule parameters:\n";
  str += "    c_inf = " + std::to_string(c_inf);
  str += ", Gamma = " + std::to_string(c_Gamma);
  str += "\n   Transition rule parameters:\n";
  str += "    B = " + std::to_string(d_B);
  str += ", pcb = " + std::to_string(d_pcb);
  str += "\n   Shear rule parameters:\n";
  str += "    mu_0 = " + std::to_string(mu_0);
  str += ", mu_1 = " + std::to_string(mu_1);
  str += "\n   Bulk rule parameters:\n";
  str += "    p_0 = " + std::to_string(K_p0);
  str += ", kappa = " + std::to_string(K_kappa);
  str += "\n   Power law exponent:\n";
  str += "    n = " + std::to_string(pl_n);
  str += "\n   Compaction function parameters:\n";
  str += "    g0 = " + std::to_string(cf_g0);
  str += ", pcinf = " + std::to_string(cf_pcinf);
  
  str += "\n";
}




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


























