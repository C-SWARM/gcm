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

#ifndef KMS_IJSS2017_templates_hpp
#define KMS_IJSS2017_templates_hpp


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






// class KMS_IJSS2017
// ******************

// Constructors
// ------------

template <int dim>
KMS_IJSS2017<dim>::KMS_IJSS2017(KMS_IJSS2017_Parameters* P )
//! Constructor with material parameters
{
  
  // Assignments
  Parameters = P;
  
  // Checks
  Parameters->Checks( true );
  
}



// Protected Methods
// -----------------


template <int dim>
double KMS_IJSS2017<dim>::a(double pc)
//! Yield function ellipses semi-axis as a function of pc
{
  return ( pc + c(pc) ) / ( 1.0 + this->Parameters->yf_alpha );
}

template <int dim>
double KMS_IJSS2017<dim>::b(double pi, double pc)
//! Yield function ellipses semi-axis as a function of pc and of pressure, pi
{
  if (pi < pi_m(pc) ) return 1.0;
  return Parameters->yf_alpha;
}

template <int dim>
double KMS_IJSS2017<dim>::betaC(double pc)
//!< Compaction function
{
  //!< Tolerance in this class
  double MODELTOL=1E-12;
  if ( fabs( pc - Parameters->cf_pcinf ) / Parameters->cf_pcinf  < MODELTOL )
    return 0;
  
  return ( Parameters->cf_g0 * ( 1.0 - pc / Parameters->cf_pcinf ) );
}

template <int dim>
double KMS_IJSS2017<dim>::betaD(double pc)
//!< Shear-induced dilatancy function
{
  return 0;
}

template <int dim>
double KMS_IJSS2017<dim>::c(double pc)
//! Cohesion as a function of pc
{
  if ( Parameters->smMbrackets )
  {
    double smf=1;
    double smM = ( 1.0 + 2.0 * atan( smf * ( pc - Parameters->d_pcb ) ) / M_PI ) * 0.5 * ( pc - Parameters->d_pcb )
               + ( 1.0 + 2.0 * atan( - smf *  Parameters->d_pcb ) / M_PI ) * 0.5 * Parameters->d_pcb ;
    
    return Parameters->c_inf * ( 1.0 - exp(- Parameters->c_Gamma *  smM ) ) ;
  }
  
  if (  pc > Parameters->d_pcb )
    return Parameters->c_inf * ( 1.0 - exp(- Parameters->c_Gamma *  ( pc - Parameters->d_pcb ) ) ) ;
  
  return 0;
}



template <int dim>
double KMS_IJSS2017<dim>::DcDpc(double pc)
//! Derivative of the Cohesion as a function of pc
{
  if ( Parameters->smMbrackets )
  {
    double smf=1;
    double tmp = smf * ( pc - Parameters->d_pcb );
    double dsmM = 0.5 + ( tmp / ( 1 + tmp * tmp ) + atan( tmp ) ) / M_PI;
    double smM = ( 1.0 + 2.0 * atan( tmp ) / M_PI ) * 0.5 * ( pc - Parameters->d_pcb )
    + ( 1.0 + 2.0 * atan( - smf *  Parameters->d_pcb ) / M_PI ) * 0.5 * Parameters->d_pcb ;
    
    return Parameters->c_inf * Parameters->c_Gamma * exp(- Parameters->c_Gamma *  smM ) * dsmM  ;
  }
  
  if (  pc > Parameters->d_pcb )
    return Parameters->c_inf * exp( - Parameters->c_Gamma * ( pc - Parameters->d_pcb ) ) * Parameters->c_Gamma;
  
  return 0;
}




template <int dim>
double KMS_IJSS2017<dim>::d(double pc)
//! Transition function as a function of pc
{
  if ( Parameters->smMbrackets )
  {
    double smf=1;
    double smM = ( 1.0 + 2.0 * atan( smf * ( pc - Parameters->d_pcb ) ) / M_PI ) * 0.5 * ( pc - Parameters->d_pcb )
                 + ( 1.0 + 2.0 * atan( - smf *  Parameters->d_pcb ) / M_PI ) * 0.5 * Parameters->d_pcb ;
    
    return 1.0 + Parameters->d_B *  smM ;
  }
  
  if (  pc > Parameters->d_pcb )
    return 1.0 + Parameters->d_B * ( pc - Parameters->d_pcb ) ;
  
  return 1.0;
}




template <int dim>
double KMS_IJSS2017<dim>::DdDpc(double pc)
//! Derivative of the transition function wrt pc
{
  if ( Parameters->smMbrackets )
  {
    double smf=1;
    double tmp = smf * ( pc - Parameters->d_pcb );
    double dsmM = 0.5 + ( tmp / ( 1 + tmp * tmp ) + atan( tmp ) ) / M_PI;
    
    return Parameters->d_B * dsmM ;
  }
  
  if ( pc > Parameters->d_pcb )
    return Parameters->d_B ;
  
  return 0.0;
}




template <int dim>
double KMS_IJSS2017<dim>::gammadot_d(double tau, double pc)
//! Equivalent shear plastic strain rate as a function of tau and pc
{
  return Parameters->flr_gamma0 * (1.0 - 1.0 / d(pc) ) * pow( tau / g_tau(pc), 1.0 / Parameters->flr_m);
}

template <int dim>
double KMS_IJSS2017<dim>::gammadot_v(double pi, double pc)
//! Equivalent volumetric plastic strain rate as a function of pi and pc
{
  if ( pi > pi_m(pc) )
    return Parameters->flr_gamma0 * pow( ( pi - pi_m(pc) ) / g_pi( pi, pc ), 1.0 / Parameters->flr_m);
  return 0;
}

template <int dim>
double KMS_IJSS2017<dim>::g_pi(double pi, double pc)
//! "Yield-like" pressure point as a function of pc
{
  return a(pc) * b(pi, pc);
}

template <int dim>
double KMS_IJSS2017<dim>::g_tau(double pc)
//! "Yield-like" "Yield-like" shear point
{
  return sqrt( 1.5 ) * a(pc) * Parameters->yf_M;
}


template <int dim>
double KMS_IJSS2017<dim>::pi_m(double pc)
//! Yield function ellipses centroid as a function of pc
{
  return a(pc) - c(pc);
}


template <int dim>
double KMS_IJSS2017<dim>::shearmodulus(double pc)
//! Shear modulus mu as a function of pc
{
  return Parameters->mu_0 + c(pc) * ( d(pc) - 1.0/d(pc) ) * Parameters->mu_1 ;
}


template <int dim>
double KMS_IJSS2017<dim>::DmuDpc(double pc)
//! Derivative of the shear modulus wrt pc
{

  double dpc = d(pc);
  return ( dpc - 1.0/dpc ) * Parameters->mu_1 * DcDpc( pc ) + c(pc) * Parameters->mu_1 * ( 1 + 1.0 / ( dpc * dpc ) ) * DdDpc( pc ) ;

}



template <int dim>
double KMS_IJSS2017<dim>::bulkmodulus(double pc)
//! Bulk modulus K as a function of pc
{
  return ( Parameters->K_p0 + c(pc) ) * ( d(pc) - 1.0/d(pc) ) / Parameters->K_kappa ;
}


template <int dim>
double KMS_IJSS2017<dim>::DbulkDpc(double pc)
//! Derivative of the bulk modulus wrt pc
{
  
  double dpc = d(pc);
  return ( dpc - 1.0/dpc ) / Parameters->K_kappa * DcDpc( pc ) + ( Parameters->K_p0 + c(pc) ) / Parameters->K_kappa * ( 1 + 1.0 / ( dpc * dpc ) ) * DdDpc( pc ) ;
  
}


template <int dim>
double KMS_IJSS2017<dim>::DalphaDpc(double pc)
//! Derivative of the coefficient alpha wrt pc
//! alpha is defined as ( 1 - c/cinf)^n / kappa
//! and is part of the definition of U(Je, pc).
{
  
  return - Parameters->pl_n * pow( 1.0 - c(pc) / Parameters->c_inf , Parameters->pl_n - 1 ) * DcDpc( pc ) / ( Parameters->K_kappa * Parameters->c_inf ) ;
  
}



template <int dim>
double KMS_IJSS2017<dim>::Dg_tauDpc(double pc)
//! Derivative of the coefficient g_tau wrt pc
//! see page 26 of the consitent tangent notes.
{
  
  return sqrt( 1.5 ) * Parameters->yf_M  / ( 1.0 + Parameters->yf_alpha )  * ( 1 + DcDpc( pc ) ) ;
  
}


template <int dim>
double KMS_IJSS2017<dim>::Dpi_mDpc(double pc)
//! Derivative of the coefficient pi_m wrt pc
//! see page 27 of the consitent tangent notes.
{
  
  return ( 1.0 - Parameters->yf_alpha * DcDpc( pc ) )   / ( 1.0 + Parameters->yf_alpha ) ;
  
}


template <int dim>
double KMS_IJSS2017<dim>::Dg_piDpc( double pi, double pc )
//! Derivative of the coefficient g_pi wrt pc
//! see page 27 of the consitent tangent notes.
{
  
  return ( 1.0 + DcDpc( pc ) ) * b( pi, pc )  / ( 1.0 + Parameters->yf_alpha );
  
}


// Public Methods
// --------------



template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017<dim>::SecondPKTensor( const FTensors& eF, const double pc )
//! Estimate of the second Piola-Kirchoff stress from a given eF
//! The volumetric part is peculiar of this model, whereas the isochoric contribution is the usual Neo-Hookean
//! with material parameters that are functions of the internal variable pc
{
  
  // tensors definitions

  using namespace ttlindexes;
  
  // FTensors Id =  ttl::Delta<2,dim,double>(1);
  FTensors Id2 = ttl::identity(i,j); // ttl::Delta<2,dim,double>(1);

  // Je at step n+1
  double eJ = det( eF ) ;
  
  // Left Cauchy-Green tensor and its inverse
  FTensors Ce;
  Ce(i,j) = eF(l,i) * eF(l,j);
  
  FTensors InvCe = ttl::inverse( Ce );
  // Inv( Ce, InvCe );
  
  // Isochoric contribution
  FTensors SPK;
  SPK(i,j) = shearmodulus(pc) * pow( eJ, -2.0/(double) dim ) * ( Id2(i,j) - Trace(Ce) / 3.0 * InvCe(i,j) );
  
  // Volumetric contribution
  double exponent = - pow( ( 1.0 - c(pc)/Parameters->c_inf  ), Parameters->pl_n  ) / Parameters->K_kappa ;
  SPK(i,j) += (
               0.5 * bulkmodulus(pc) * ( eJ * log(eJ) + eJ - 1.0 )
               + ( c(pc) - ( Parameters->K_p0 + c(pc) ) * pow( eJ, exponent ) )
               )
  * InvCe(i,j) ;
  
  return SPK;
  
}



template <int dim>
double KMS_IJSS2017<dim>::HardeningLaw(double pc)
// Hardening law, returns Log(jp) as a function of pc
{
  double a1=this->Parameters->hr_a1, a2=this->Parameters->hr_a2, l1=this->Parameters->hr_Lambda1, l2=this->Parameters->hr_Lambda2 ;
  return - a1 * exp( - l1 / pc ) - a2 * exp( - l2 / pc ) ;
}


template <int dim>
double KMS_IJSS2017<dim>::DHardeningLawDpc(double pc)
// Derivative of the Hardening law as a function of pc
{
  double a1=this->Parameters->hr_a1, a2=this->Parameters->hr_a2, l1=this->Parameters->hr_Lambda1, l2=this->Parameters->hr_Lambda2 ;
  return - ( a1 * l1 * exp( - l1/ pc ) + a2 * l2 * exp( - l2/ pc ) ) / ( pc * pc ) ;
}




template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017<dim>::KirchoffStressTensor( const FTensors& Fp, const FTensors& Fe, const FTensors& S )
//! Estimate of the Kirchoff stress from a given pF, eF, S
{
  using namespace ttlindexes;
  return det(Fp) * Fe(i,j) * S(j,k) * Fe(l,k);
}


template <int dim>
ttl::Tensor<2, dim, double> KMS_IJSS2017<dim>::CauchyStressTensor( const FTensors& Fe, const FTensors& S )
//! Estimate of the Cauchy stress from a given eF, S
{
  using namespace ttlindexes;
  return Fe(i,j) * S(j,k) * Fe(l,k) / det(Fe) ;
}







#endif /* KMS_IJSS2017_templates_hpp */
