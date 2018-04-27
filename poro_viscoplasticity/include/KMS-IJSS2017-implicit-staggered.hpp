//
//  KMS-IJSS2017-implicit-staggered.hpp
//  ttl-learning
//
//  Created by alberto salvadori on 2/1/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

/*!
 
 \file KMS-IJSS2017.h
 \brief Template code file for the KMS-IJSS2017 model
 
 This file contains template classes and methods for the staggered implicit integration of the poroviscoplastic model desccribed in the paper
 "A. Krairi, K. Matouš, A. Salvadori, A poro-viscoplastic constitutive model for granular materials at finite strain, submitted to IJSS (2017)".
 
 */


#ifndef KMS_IJSS2017_implicit_staggered_h
#define KMS_IJSS2017_implicit_staggered_h



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


// class KMS_IJSS2017_Implicit_BE_Staggered
// staggered backward euler integrator
// ****************************************



template <int dim>
KMS_IJSS2017_Implicit_BE_Staggered<dim>::KMS_IJSS2017_Implicit_BE_Staggered(KMS_IJSS2017_Parameters* P, KMS_IJSS2017_IO* InpOut, TimeIntegrationDataManager* TIDM )
: KMS_IJSS2017_Implicit_BE<dim>( P, InpOut, TIDM )
//! Constructor with material parameters and IO
{}


// IO - methods
// ------------

template <int dim>
void KMS_IJSS2017_Implicit_BE_Staggered<dim>::AsAString( std::string& str, bool Verbose )
//! This method prints the model features as a string
{
  
  if ( Verbose )
  {
    // Parameters
    if (this->Parameters == NULL )
      str += " No parameters have been set, yet. \n";
    else
      this->Parameters->AsAString(str);
    
    // Integrators
    str += "\n";
    str += "  KMS_IJSS2017 Staggered Implicit Integrator for dimension ";
    str += std::to_string(dim);
    str += ". \n";
  }
  
}


template <int dim>
unsigned KMS_IJSS2017_Implicit_BE_Staggered<dim>::explicit_integrator( const FTensors& updF, const double dt, bool Verbose )
//! Updates stresses and internal variables, given Fnp1
/*!
 // Explicit integration algorithm
 // 1 - Assign DeltaF
 // 2 - Compute strain rates
 // 3 - Compute plasting stretching Dp
 // 4 - Compute pFnp1
 // 5 - Compute pcnp1
 // 6 - Compute eFnp1
 // 7 - Compute Snp1
 // 8 - Compute KSnp1
 // 9 - Compute sigmanp1
 
 */
{
  
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;
  
  // Kirchoff stress definition
  ttl::Tensor<2, dim, double> Id2 = ttl::identity( i,j ); // ttl::Delta<2,dim,double>(1);
  
  FTensors Kirchoffn;
  FTensors DevKirchoffn;
  
  // Kirchoff stress construction
  double Jn = det( this->Fn ), eJn = det( this->eFn );
  Kirchoffn(i,j) = Jn/eJn * this->eFn(i,k) * this->Sn(k,l) * this->eFn(j,l);
  
  // spherical part -
  // note that p is taken as positive in compression, whence the sign - before the Trace operator
  double pi = - Trace( Kirchoffn ) / 3.0 ;
  //deviatoric part
  DevKirchoffn = dev( Kirchoffn );
  //Dev( Kirchoffn, DevKirchoffn );
  // Frobenius norm
  double FrobNormDevKirchoffn = FrobeniusNorm( DevKirchoffn );
  
  // 1
  this->Fnp1 = updF;
    
  // 2
  double gmdv = this->gammadot_v(pi, this->pcn);
  double gmdd = this->gammadot_d(FrobNormDevKirchoffn, this->pcn);
  
  if ( std::fpclassify(gmdv) == FP_ZERO && std::fpclassify(gmdd) == FP_ZERO)
  // if ( fabs(gmdv) < this->ZERO && fabs(gmdd) < this->ZERO )
  {    
    this->Dpn = Id2(i,j) * 0 ; // ttl::Delta<2,dim,double>(0);
    this->pFnp1(i,j) =  this->pFn(i,j);
    this->pcnp1 = this->pcn;
    
  }
  else
  {    
    // 3
    FTensors Nn = dev( this->Sn );
    double NnNorm= FrobeniusNorm( Nn );
    
    if ( NnNorm > this->INTEGRATOR_TOL )
      this->Dpn(i,j) = gmdd / NnNorm * Nn(i,j) + 1.0/3.0 * ( this->betaD( this->pcn ) * gmdd - this->betaC( this->pcn ) * gmdv ) * Id2(i,j);
    else
      this->Dpn(i,j) = 1.0/3.0 * ( this->betaD( this->pcn ) * gmdd - this->betaC( this->pcn ) * gmdv ) * Id2(i,j);
        
    // 4
    this->pFnp1(i,j) =  this->pFn(i,j) + this->Dpn(i,l) * this->pFn(l,j) * this->TimeIntegrationData->Deltat();
    
    // 5
    this->FindpcFromJpAtStepnP1( Verbose );
  }

  // 6
  FTensors buffer = ttl::inverse( this->pFnp1 );
  // Inv( this->pFnp1, buffer );
  this->eFnp1(i,j) =  this->Fnp1(i,l) * buffer(l,j);
  
  // 7
  this->SecondPKTensorAtStepnP1( Verbose );
  
  // 8
  this->KSnp1 = this->KirchoffStressTensor( this->pFnp1, this->eFnp1, this->Snp1 );
  
  // 9
  this->sigmanp1 = this->CauchyStressTensor( this->eFnp1, this->Snp1 );
  
  return 0;
  
}


template <int dim>
unsigned KMS_IJSS2017_Implicit_BE_Staggered<dim>::FindpcFromJp( const double logJp, double& pcr, double PCTOL, bool Verbose )
//! NR scheme to estimate pc at step (n+1) from Jp at the same step
//! Returns the number of iterations required to convergence, whereas the new outcome is stored in pcr and returned as such
//! The input data is contained in logJpr that must be updated before this function is called
{  
  if(logJp>0)
    return 0;
  
  // if pcr was set to pcinf, there is nothing to do since pc is monotonically increasing.

  if ( this->pcEQpc_inf )
    return 0;
  
  // otherwise look for pc
  
  double pcr0 = pcr;
  double NRTOL = 1e-6;
  
  // Initial checks
  
  double a1=this->Parameters->hr_a1, a2=this->Parameters->hr_a2, l1=this->Parameters->hr_Lambda1, l2=this->Parameters->hr_Lambda2;
  
  double fpcr = this->HardeningLaw( pcr ) - logJp ;
  
  if ( std::fpclassify( fpcr ) == FP_ZERO )
    return 0;
  
  if ( Verbose && ( fpcr < 0 ) )
  {
    std::cout << " WARNING: KMS_IJSS2017_Implicit_BE_Staggered<dim>::FindpcFromJp( ) - Unexpected sign for f(pc(r)): it is negative, which means that pc decreased.";
    std::cerr << " Code did not abort but outcomes might be wrong. \n";
  }
  
  double DfDpcr = ( a1 * l1 * exp( - l1/ this->pcn ) + a2 * l2 * exp( - l2/ this->pcn ) ) / ( this->pcn * this->pcn );
  
  // Algorithm selection, either NR or bisections
  
  unsigned it=0;
  double computer_zero = 1.0e-15;
  
  if(fabs( DfDpcr ) < NRTOL )  
    // the first derivative is almost zero, and the NR scheme becomes very bad conditioned.
    // the zero is found via bisection method
  {
    if ( Verbose )
    {
      std::cout << " WARNING: KMS_IJSS2017_Implicit_BE_Staggered<dim>::FindpcFromJp( ) - The first derivative is almost zero at pc(n). pc(n+1) will be estimated via bisection method. ";
    }
    
    double a = pcr0, b = a, fpcrp1 = 1, c=a;
    
    // seeking for the upper extreme of the bounding interval
    if ( fpcr > 0 )
    {  
      do
      {
        b *= 2;
        fpcrp1 = this->HardeningLaw( b ) - logJp ;
      } while ( fpcrp1 > 0 );
    }  
    else
    {  
      do
      {
        a /= 2;
        double H = this->HardeningLaw( a );
        fpcrp1 = H - logJp ;
        
 //       if(fabs(H)<1.0e-30)
 //         break;
        
      } while ( fpcrp1 < 0 );
    }
    
    
    // at this point we have fpcr>0 and fpcrp1 <0
    // the bisection method can run
    unsigned maxIt=500;
    double fpccheck = 1, fa=1;
    
    while (it < maxIt)
    {
      it++;
      c = 0.5 * (a+b);
      
      fa = this->HardeningLaw( a ) - logJp ;
      fpccheck = this->HardeningLaw( c ) - logJp ;
      
      if ( fabs( fpccheck ) < computer_zero )
        break;
      else if ( fa * fpccheck > 0 )
        a = c;
      else
        b = c;
    }
    
    pcr = c;
    
    if ( it == maxIt )
    {
      //this->IO->log() << " WARNING: KMS_IJSS2017_Implicit_BE_Staggered<dim>::FindpcFromJp( ) - Convergence was not achieved in the bisection scheme after " << maxIt << " iterations." << std::flush;
      //this->IO->log() << " Code did not abort but outcomes might be wrong. \n" << std::flush;
    }
  }
  

    // the NR method is used to find pcr
  {
    unsigned maxIt=100;
//    pcr = pcr0;
//    pcr = (l1+l2)/2.0;
    
    double R = this->HardeningLaw( pcr ) - logJp;
    double norm   = sqrt(R*R);
    double norm_0 = norm;
    
    if(norm_0<computer_zero)
      norm_0 = computer_zero;
    
    while (it < maxIt) {
      
      it++;
      double dpc = ( this->HardeningLaw( pcr ) - logJp ) * ( pcr * pcr );
      dpc /= a1 * l1 * exp( - l1 / pcr ) + a2 * l2 * exp( - l2 / pcr ) ;
      pcr += dpc;
      
      if ( fabs(  this->HardeningLaw( pcr ) - logJp  ) < computer_zero)
        break;
    }
    
    if ( it == maxIt )
    {
      //this->IO->log() << " WARNING: KMS_IJSS2017_Implicit_BE_Staggered<dim>::FindpcFromJp( ) - Convergence was not achieved in the Newton-Raphson scheme after " << maxIt << " iterations." << std::flush;
      //this->IO->log() << " Code did not abort but outcomes might be wrong. \n" << std::flush;
    }
  }
    
  if(pcr<pcr0)
    pcr = pcr0;

  return it;
}

template <int dim>
unsigned KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose )
//! Updates stresses and internal variables, given Fnp1
//! The algorithm loops until tolerance STAGGEREDTOL is achieved.
//! The initial data are all those at step n, and Fnp1 = updF
//! The algorithm evaluates pFnp1, eFnp1, pcnp1
//!
//! Staggered imlicit integration algorithm
//! 1 - Initialize: Assign Fnp1, deltat
//! 2 - Solve for r until convergence
//! 3 - Compute pFnp1
//! 4 - Compute eFnp1
//! 5 - Compute Snp1
//! 6 - Compute KSnp1
//! 7 - Compute sigmanp1
//!
//! All prints to std::cout have been removed, the function VerboseStepUpdate( const FTensors& updF, const double dt, bool Verbose )
//! shall be used in debug mode
{
  using namespace ttlindexes;  
  FTensors dF = this->Fnp1(i,j) - this->Fnp1(i,j);
    
  if(FrobeniusNorm(dF)<1.0e-12 || dt<1.0e-12)
  {
    return this->explicit_integrator(updF, dt, false);
  }  
  
  // 1: r = 0
  
  // Updated tensors and their clean up.
  // The need of cleaning tensors stands in the fact that very small off-diagonal terms (of the order E-90)
  // cause issues in the Mechanical Model Integration. Rather, if they are set to zero, no issues whatsoever
  // arise. The check can be done on F since the diagonal term should be of magnitude 1.

  this->Fnp1 = updF;
  this->CleanFTensors( this->Fnp1 );
  
  FTensors InvFnp1 = ttl::inverse( this->Fnp1 );
  
  // Intermediate tensors initialization
  // eFr is initialized with the amount at the end of step n
  // pFr is initialized according to the multiplicative decomposition F = eF pF
  //     i.e. pFr = inverse(eFr) * Fnp1
  // Sr  is initialized with the amount at the end of step n, since eF and pc are initialized in the same way
  
  double tol = 1.0e-12;
  
  FTensors eFr = this->eFn;
  this->CleanFTensors( eFr );

  FTensors pFr = ttl::inverse(eFr)( i,k ) * this->Fnp1( k,j );
  
  FTensors Sr = this->Sn;
  FTensors Mr;
  
  Mr( i,j ) = this->pFn( i,k ) * InvFnp1( k,l ) * eFr( l,j );
  
  double pcr = this->pcn;
  
  FTensors InvMr = {};
  
  
  // 2: Solve for r until convergence.
  
  // external loop with iterator it
  // it is broken when convegence is gained for M and pc
  
  unsigned it=0;
  while (it < MAX_N_OF_STAGGERED_ITERATIONS) {
    
    it++;
    
    // staggered solution for M
    // internal loop wit iterator itMr
    // it is broken when convegence is gained for M
    
    unsigned itMr=0;
    while (itMr < MAX_N_OF_STAGGERED_ITERATIONS) {
      
      ttl::Tensor<2, dim, double> rm = this->RM( dt, this->pFn, this->Fnp1, eFr, Sr, Mr, pcr );
      
      // Convergence check on the residual for M
      
      if ( FrobeniusNorm( rm ) / FrobeniusNorm( Mr )  < tol  )
        break;
      
      // if convergence is not achieved, iterate
      
      itMr++;
      
      ttl::Tensor<4, dim, double> drdm = this->DRMDM( dt, this->pFn, this->Fnp1, eFr, Sr, Mr, pcr );
      ttl::Tensor<2, dim, double> dm = {};
      
      try
      {
      
        dm = ttl::solve( drdm, rm );
      
      }
      catch (int i)
      {

        ttlsolveExceptionHandling( drdm, rm );
      
      }
      

      Mr( i,j ) = Mr( i,j ) - dm( i,j );

      // 3-5 : pF, eF, Sr shall be updated:
      
      // 3 - Compute pFr

      try{
        InvMr = ttl::inverse( Mr );
      }
      catch(const int ttlerr){
        throw 1;
      }      
      
      pFr( i,k ) = InvMr( i,j ) * this->pFn( j,k );
      
      // 4 - Compute eFr
      
      FTensors InvpFr = ttl::inverse( pFr );

      eFr( i,k ) = this->Fnp1( i,j ) * InvpFr( j,k );
      
      // 5 - Compute Sr
      
      Sr = this->SecondPKTensor( eFr, pcr );
      
    }
    
    // staggered solution for pc
    // the numerical solution is estimated in the method FindpcFromJp, that makes use of the Newton-Raphson
    // scheme when the problem is well conditioned, and of the bisection method when the problem is ill conditioned
    
    double logJpr = log( det( pFr ) );
      FindpcFromJp( logJpr, pcr, tol, false ); // Verbose );
    
    // convergence check for the residual of the flow rule,
    // to exit the whole staggered algorithm
    
    ttl::Tensor<2, dim, double> rm = this->RM( dt, this->pFn, this->Fnp1, eFr, Sr, Mr, pcr );
    
    if ( FrobeniusNorm( rm ) / FrobeniusNorm( Mr )  < tol  )
      break;
    
  }
  
  if ( it == MAX_N_OF_STAGGERED_ITERATIONS )
  {
    //this->IO->log() << " WARNING: KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose ) - Convergence was not achieved for it in the staggered scheme after " << MAX_N_OF_STAGGERED_ITERATIONS << " iterations.";
    //this->IO->log() << " Code did not abort but outcomes might be wrong. \n";
  }
  
  
  // if pcr is numerically very close to or even above pcinf, it is set to pcinf
  // so that the controls take care of possible overflows
  
  if ( this->pcEQpc_inf == false)
    if ( ( 1 - pcr / this->Parameters->cf_pcinf ) < 1E-6 )
    {
      pcr =  this->Parameters->cf_pcinf;
      this->pcEQpc_inf = true;

//      this->IO->log() << "\n  Time = " << std::fixed << std::setprecision(3) << this->TimeIntegrationData->CurrentTime() << ", step = " << this->TimeIntegrationData->TimeStep() ;
//      this->IO->log() << ", pcr set equal to pc_inf \n\n" << std::flush;
//      this->IO->log() << " / pcr set equal to pc_inf" << std::flush;
    }
  
  // if pc overcomes pcb, print it out.
  if (
      ( this->pcn < this->Parameters->d_pcb )
      &&
      ( pcr > this->Parameters->d_pcb )
      )
  {
//    this->IO->log() << "\n  Time = " << std::fixed << std::setprecision(3) << this->TimeIntegrationData->CurrentTime() << ", step = " << this->TimeIntegrationData->TimeStep() ;
//    this->IO->log() << ", pcr (" << std::setw(20) << std::setprecision(15) << pcr
//                    << ")overcomes pcb (" << std::setw(20) << std::setprecision(15) << this->Parameters->d_pcb << ")\n\n" << std::flush;
//    this->IO->log() << " / pcr (" << std::setw(20) << std::setprecision(15) << pcr
//    << ") overcomes pcb (" << std::setw(20) << std::setprecision(15) << this->Parameters->d_pcb << ")" << std::flush;
  }

  // The algorithm converged for M and for pc.
  // pc, pF, eF, and S at step n+1 are equal to the r-ones
  
  this->pcnp1 = pcr;
  this->pFnp1 = pFr;
  this->eFnp1 = eFr;
  this->Snp1 = Sr;
  
  //! 6 - Compute KSnp1
  
  this->KSnp1 = this->KirchoffStressTensor( this->pFnp1, this->eFnp1, this->Snp1 );
  
  //! 7 - Compute sigmanp1
  
  this->sigmanp1 = this->CauchyStressTensor( this->eFnp1, this->Snp1 );
  
  // return the number of iterations
  
  return it;
  
}





template <int dim>
unsigned KMS_IJSS2017_Implicit_BE_Staggered<dim>::VerboseStepUpdate( const FTensors& updF, const double dt, bool Verbose )
//! Updates stresses and internal variables, given Fnp1
//! The algorithm loops until tolerance STAGGEREDTOL is achieved.
//! The initial data are all those at step n, and Fnp1 = updF
//! The algorithm evaluates pFnp1, eFnp1, pcnp1
//!
//! Staggered imlicit integration algorithm
//! 1 - Initialize: Assign Fnp1, deltat
//! 2 - Solve for r until convergence
//! 3 - Compute pFnp1
//! 4 - Compute eFnp1
//! 5 - Compute Snp1
//! 6 - Compute KSnp1
//! 7 - Compute sigmanp1
//!
{
  
  using namespace ttlindexes;
  
  // 1: r = 0
  
  // Updated tensors and their clean up.
  // The need of cleaning tensors stands in the fact that very small off-diagonal terms (of the order E-90)
  // cause issues in the Mechanical Model Integration. Rather, if they are set to zero, no issues whatsoever
  // arise. The check can be done on F since the diagonal term should be of magnitude 1.
  
  this->Fnp1 = updF;
  this->CleanFTensors( this->Fnp1 );
  
  FTensors InvFnp1 = ttl::inverse( this->Fnp1 );
  
  // Intermediate tensors initialization
  // eFr is initialized with the amount at the end of step n
  // pFr is initialized according to the multiplicative decomposition F = eF pF
  //     i.e. pFr = inverse(eFr) * Fnp1
  // Sr  is initialized with the amount at the end of step n, since eF and pc are initialized in the same way
  
  FTensors eFr = this->eFn;
  this->CleanFTensors( eFr );
  
  FTensors pFr = ttl::inverse(eFr)( i,k ) * this->Fnp1( k,j );

  FTensors Sr = this->Sn;
  FTensors Mr;
  Mr( i,j ) = this->pFn( i,k ) * InvFnp1( k,l ) * eFr( l,j );
  
  double pcr = this->pcn;
  
  FTensors InvMr = {};
  
  if ( Verbose )
  {
    std::string Fstr="   Initialization of KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose ) \n";
    std::cout << Fstr;
    Fstr = "";
    PrintInVectorForm( this->Fnp1 , Fstr );
    std::cout << "    F(n+1) = " << Fstr;
    Fstr = "";
    PrintInVectorForm( pFr , Fstr );
    std::cout << "    pFr    = " << Fstr;
    Fstr = "";
    PrintInVectorForm( eFr , Fstr );
    std::cout << "    eFr    = " << Fstr;
    Fstr = "";
    PrintInVectorForm( Sr , Fstr );
    std::cout << "    Sr     = " << Fstr;
    Fstr = "";
    PrintInVectorForm( Mr , Fstr );
    std::cout << "    Mr     = " << Fstr;
    Fstr = "";
    std::cout << "    pcr    = " << std::setw(20) << std::setprecision(15) << pcr ;
    std::cout << ", pcb    = " << std::setw(20) << std::setprecision(15) << this->Parameters->d_pcb << "\n";
  }
  
  
  // 2: Solve for r until convergence.
  
  if ( Verbose )
  {
    std::cout << "\n   Staggered iterations \n";
  }
  
  // external loop with iterator it
  // it is broken when convegence is gained for M and pc
  
  unsigned it=0;
  while (it < MAX_N_OF_STAGGERED_ITERATIONS) {
    
    it++;
    
    // staggered solution for M
    // internal loop wit iterator itMr
    // it is broken when convegence is gained for M
    
    unsigned itMr=0;
    while (itMr < MAX_N_OF_STAGGERED_ITERATIONS) {
      
      ttl::Tensor<2, dim, double> rm = this->RM( dt, this->pFn, this->Fnp1, eFr, Sr, Mr, pcr );

      if ( Verbose )
      {
        std::cout << "    it = " << it << ", itMr = " << itMr  << ", || rm || = " << std::setw(20) << std::setprecision(15) << FrobeniusNorm( rm )
                  << ", || Mr || = " << std::setw(20) << std::setprecision(15) << FrobeniusNorm( Mr ) << std::flush;

      }

      // Convergence check on the residual for M
      
      if ( FrobeniusNorm( rm ) / FrobeniusNorm( Mr )  < STAGGEREDTOL  )
        break;
      
      // if convergence is not achieved, iterate
      
      itMr++;
      
      // ttl::Tensor<4, dim, double> drdm = this->DRMDM( this->Deltat, this->pFn, this->Fnp1, eFr, Sr, Mr, pcr );
      ttl::Tensor<4, dim, double> drdm = this->DRMDM( dt, this->pFn, this->Fnp1, eFr, Sr, Mr, pcr );
      
      if ( Verbose )
      {
        
        std::string drdmstr="";
        PrintInVectorForm( drdm, drdmstr );
        
        std::cout << "\n" << "     drdm  = " << drdmstr ;
        std::cout << ", || drdm || = " << std::setw(20) << std::setprecision(15) << FrobeniusNorm(drdm) << std::flush;
        
        std::string rmstr="";
        PrintInVectorForm( rm, rmstr );
        
        std::cout << "\n" << "     rm  = " << rmstr ;
        
      }
      
      //ttl::Tensor<2, dim, double> dm;
      //int sysol = ttl::solve( drdm, rm, dm );
      
      ttl::Tensor<2, dim, double> dm = {};
      
      try
      {
        
        dm = ttl::solve( drdm, rm );
        
      }
      catch (int i)
      {
        
        ttlsolveExceptionHandling( drdm, rm );
        
      }


      if ( Verbose )
      {
        std::string Fstr = "";
        PrintInVectorForm( dm , Fstr );
        std::cout << "\n     dm     = " << Fstr;
        Fstr = "";

        std::cout << ", || dm || = " << std::setw(20) << std::setprecision(15) << FrobeniusNorm(dm) << "\n";
        // std::cout << "sysol = " << sysol << "\n";
        
      }

      ttl::Tensor<2, dim, double> solvercheck = drdm( i,j,k,l ) * dm( k,l ) - rm ( i,j );
      
      if ( Verbose )
      {
        std::string Fstr = "";
        PrintInVectorForm( solvercheck , Fstr );
        std::cout << "\n     solvercheck     = " << Fstr;
        Fstr = "";
        
        std::cout << ", || solvercheck || = " << std::setw(20) << std::setprecision(15) << FrobeniusNorm(solvercheck) << "\n";
        
      }

      Mr( i,j ) = Mr( i,j ) - dm( i,j );
      
      if ( Verbose )
      {
        
        std::string Fstr = "";
        PrintInVectorForm( Mr , Fstr );
        std::cout << "             Mr     = " << Fstr;
        Fstr = "";
      }
      
      // 3-5 : pF, eF, Sr shall be updated:

      // 3 - Compute pFr
      
      InvMr = ttl::inverse( Mr );
      pFr( i,k ) = InvMr( i,j ) * this->pFn( j,k );
      
      // 4 - Compute eFr
      
      FTensors InvpFr = ttl::inverse( pFr );
      eFr( i,k ) = this->Fnp1( i,j ) * InvpFr( j,k );
      
      // 5 - Compute Sr
      
      Sr = this->SecondPKTensor( eFr, pcr );

    }

    if ( Verbose && itMr == MAX_N_OF_STAGGERED_ITERATIONS )
    {
      std::cerr << " WARNING: KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose ) - Convergence in itMr was not achieved in the staggered scheme after " << MAX_N_OF_STAGGERED_ITERATIONS << " iterations.";
      std::cerr << " Code did not abort but outcomes might be wrong. \n";
    }

    // After Mr converged, the related tensors shall be updated.
    // They are pF, eF, Sr
    // other tensors are not necessary to be updated at this time and they will not.
    
    if ( Verbose )
    {
      std::string Fstr="\n\n    Staggered solution for Mr, pFr, eFr, Sr. \n";
      std::cout << Fstr;
      Fstr = "";
      PrintInVectorForm( Mr , Fstr );
      std::cout << "     Mr     = " << Fstr;
      Fstr = "";
      PrintInVectorForm( pFr , Fstr );
      std::cout << "     pFr    = " << Fstr;
      Fstr = "";
      PrintInVectorForm( eFr , Fstr );
      std::cout << "     eFr    = " << Fstr;
      Fstr = "";
      PrintInVectorForm( Sr , Fstr );
      std::cout << "     Sr     = " << Fstr;
      Fstr = "";
      std::cout << "     pcr    = " << std::setw(20) << std::setprecision(15) << pcr ;
      std::cout << ", pcb    = " << std::setw(20) << std::setprecision(15) << this->Parameters->d_pcb << "\n";
    }

    // staggered solution for pc
    // the numerical solution is estimated in the method FindpcFromJp, that makes use of the Newton-Raphson
    // scheme when the problem is well conditioned, and of the bisection method when the problem is ill conditioned
    
    double logJpr = log( det( pFr ) );
    unsigned convitforpc = FindpcFromJp( logJpr, pcr, STAGGEREDTOL, false ); // Verbose );
    
    if ( Verbose )
    {
      std::cout << "\n    it = " << it << ", convergence achieved for pc in  = " << convitforpc  << " iterations. pcr = " << std::setw(20) << std::setprecision(15) << pcr;
      std::cout << ", || Rpc || = " << std::setw(20) << std::setprecision(15) << fabs( this->HardeningLaw( pcr ) - logJpr ) << "\n";
    }

    // convergence check for the residual of the flow rule,
    // to exit the whole staggered algorithm
    
    ttl::Tensor<2, dim, double> rm = this->RM( dt, this->pFn, this->Fnp1, eFr, Sr, Mr, pcr );
    
    if ( Verbose )
    {
      std::cout << "    it = " << it << ", || rm || = " << std::setw(20) << std::setprecision(15) << FrobeniusNorm(rm) ;
      std::cout << ", || Mr || = " << std::setw(20) << std::setprecision(15) << FrobeniusNorm(Mr) << "\n\n";
    }
    
    if ( FrobeniusNorm( rm ) / FrobeniusNorm( Mr )  < STAGGEREDTOL  )
      break;
  
  }
  
  if ( Verbose && it == MAX_N_OF_STAGGERED_ITERATIONS )
  {
    std::cerr << " WARNING: KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose ) - Convergence was not achieved for it in the staggered scheme after " << MAX_N_OF_STAGGERED_ITERATIONS << " iterations.";
    std::cerr << " Code did not abort but outcomes might be wrong. \n";
  }

    // if pcr is numerically very close to or even above pcinf, it is set to pcinf
    // so that the controls take care of possible overflows
  
  if ( this->pcEQpc_inf == false)
    if ( ( 1 - pcr / this->Parameters->cf_pcinf ) < 1E-6 )
    {
      pcr =  this->Parameters->cf_pcinf;
      this->pcEQpc_inf = true;
      
        //      this->IO->log() << "\n  Time = " << std::fixed << std::setprecision(3) << this->TimeIntegrationData->CurrentTime() << ", step = " << this->TimeIntegrationData->TimeStep() ;
        //      this->IO->log() << ", pcr set equal to pc_inf \n\n" << std::flush;
      std::cout << " / pcr set equal to pc_inf" << std::flush;
    }
  
    // if pc overcomes pcb, print it out.
  if (
      ( this->pcn < this->Parameters->d_pcb )
      &&
      ( pcr > this->Parameters->d_pcb )
      )
  {
      //    this->IO->log() << "\n  Time = " << std::fixed << std::setprecision(3) << this->TimeIntegrationData->CurrentTime() << ", step = " << this->TimeIntegrationData->TimeStep() ;
      //    this->IO->log() << ", pcr (" << std::setw(20) << std::setprecision(15) << pcr
      //                    << ")overcomes pcb (" << std::setw(20) << std::setprecision(15) << this->Parameters->d_pcb << ")\n\n" << std::flush;
    std::cout << " / pcr (" << std::setw(20) << std::setprecision(15) << pcr
    << ") overcomes pcb (" << std::setw(20) << std::setprecision(15) << this->Parameters->d_pcb << ")" << std::flush;
  }
  

  
  if ( Verbose )
  {
    std::cout << "\n    Convergence achieved for pc and M  \n";
  }

  // The algorithm converged for M and for pc.
  // pc, pF, eF, and S at step n+1 are equal to the r-ones
  
  this->pcnp1 = pcr;
  this->pFnp1 = pFr;
  this->eFnp1 = eFr;
  this->Snp1 = Sr;
  
  if ( Verbose )
  {
    std::cout << "    pc = " << std::setw(20) << std::setprecision(15) << this->pcnp1 << "\n";

    std::string Fstr="";
    PrintInVectorForm( this->Snp1 , Fstr );
    std::cout << "    Second Piola-Kirchoff S(n+1)=" << Fstr ;
  }
   

  //! 6 - Compute KSnp1
  
  this->KSnp1 = this->KirchoffStressTensor( this->pFnp1, this->eFnp1, this->Snp1 );
  
  if ( Verbose )
  {
    std::string Fstr="";
    PrintInVectorForm( this->KSnp1 , Fstr );
    std::cout << "    Kirchoff stress S(n+1)=" << Fstr;
  }

  //! 7 - Compute sigmanp1
  
  this->sigmanp1 = this->CauchyStressTensor( this->eFnp1, this->Snp1 );
  
  if ( Verbose )
  {
    std::string Fstr="";
    PrintInVectorForm( this->sigmanp1 , Fstr );
    std::cout << "    Cauchy stress S(n+1)=" << Fstr ;
  }
  
  // return the number of iterations
  
  return it;

}






template <int dim>
void KMS_IJSS2017_Implicit_BE_Staggered<dim>::ttlsolveExceptionHandling( const ttl::Tensor<4, dim, double>& ttlsysmat, const FTensors& ttlsysvec)
//! Updates stresses and internal variables, given Fnp1
//! The algorithm loops until tolerance STAGGEREDTOL is achieved.
//! The initial data are all those at step n, and Fnp1 = updF
//! The algorithm evaluates pFnp1, eFnp1, pcnp1
//!
//! Staggered imlicit integration algorithm
//! 1 - Initialize: Assign Fnp1, deltat
//! 2 - Solve for r until convergence
//! 3 - Compute pFnp1
//! 4 - Compute eFnp1
//! 5 - Compute Snp1
//! 6 - Compute KSnp1
//! 7 - Compute sigmanp1
//!
{
  
  //this->IO->log() << "\n\n  ----------------------------- \n";
  //this->IO->log() << "  ttl::solve Exception Handling \n";
  //this->IO->log() << "  the code is aborting at \n\n  ";
  //this->IO->log() << "  timestep = " <<  this->TimeIntegrationData->TimeStep()  << ", time = " <<  this->TimeIntegrationData->CurrentTime() ;
  //this->IO->log() << "\n\n  because lapack cannot solve the linear system \n";
  //this->IO->log() << "  within the step update in the Newton-Raphson scheme of the  \n";
  //this->IO->log() << "  KMS_IJSS2017_Implicit_BE_Staggered integrator. Here are the system matrix \n";
  //this->IO->log() << "  and the right hand side before quitting. \n";
  
//  typedef mtl::mat::dense2D<double> Matrix;
//  typedef mtl::dense_vector<double> Vector;
//
//  // Matrix analysis - using mtl
//  Matrix SysMat( dim*dim, dim*dim );
//  
//  unsigned r=0;
//  for (unsigned i=0; i<dim*dim; i++)
//    for (unsigned j=0; j<dim*dim; j++)
//      SysMat[i][j] = ttlsysmat.data[r++];
//
//  this->IO->log() << "\n  system matrix: \n" << SysMat << "\n";
//  
//  // Vector analysis - using mtl
//  Vector SysVec( dim*dim );
//  
//  r=0;
//  for (unsigned i=0; i<dim*dim; i++)
//      SysVec[i] = ttlsysvec.data[r++];
//  
//  this->IO->log() << "  right hand side:   " << SysVec << "\n\n\n";
//  
//  std::exit(-1);
  
}



template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE_Staggered<dim>::IterativeDMDF( ttl::Tensor<4, dim, double> initdmdf,
                                                                                   const double dt,
                                                                                   const bool Verbose
                                                                                   )
//! This function evaluates the fourth order tensor DM/DF for the KMS_IJSS2017 model
//! in an iterative way. It is used for debug only.
//! This method expects that the Step has already been updated elsewhere, since it uses quantities at step n+1
//! that MUST have been updated before
//!
//! 1 - Calculate intermediate operators
//! 2 - Calculate some derivatives wrt pc
//! 3 - Calculate blocks 2.1, 2.2, 2.3, 2.4 and right hand sides
//! 4 - Solve the linear system and return
{
  
  
  using namespace ttlindexes;
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'y'> y;
  
  double onethird = 1.0 / 3.0;
  
  
  // 1 - Calculate intermediate operators
  
  // tensor M = pFn * InvpFnp1
  FTensors InvpFnp1 = ttl::inverse( this->pFnp1 );
  FTensors M = this->pFn( i,k ) * InvpFnp1( k,j );
  
  // Zwx has been defined at page 17b of the consistent tangent operator notes
  FTensors Zwx = this->UpdateZwx( M );
  
  // Lambdazywx has been defined at page 16 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> Lambdazywx;
  ttl::Tensor<4, dim, double> dsdm = this->DSDM( this->Fnp1, this->eFnp1, this->pFnp1, this->pcnp1 );
  Lambdazywx( z,y,w,x ) = dsdm( z,y,w,x ) + this->DSDpc() ( z,y ) * Zwx( w,x );
  
  // Vzykl has been defined at page 16 of the consistent tangent operator notes
  // Note that dSdF is the partial derivative and differs from DSDF
  ttl::Tensor<4, dim, double> Vzykl = this->dSdF( this->Fnp1, this->eFnp1, this->pFnp1, this->pcnp1 );
  
  // Minorcdkl and Majorcdwx have been defined at pages 16-17 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> Minorcdkl = this->Minor_cdkl( this->Fnp1, this->eFnp1, this->pFnp1, this->Snp1, M, Vzykl );
  ttl::Tensor<4, dim, double> Majorcdwx = this->Major_cdwx( this->Fnp1, this->eFnp1, this->pFnp1, this->Snp1, M, Lambdazywx );
  
  // Phiabzy has been defined at page 18 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> Phiabzy = this->Phi_abzy( this->Snp1, false );
  
  
  // 2 - Calculate some derivatives wrt pc
  
  // taubar is the Frobenius norm of the deviatoric part of the Kirchoff stress.
  // Note that this->KSnp1 has been updated by the method
  // KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose )
  FTensors tau = dev( this->KSnp1 );
  double taubar = FrobeniusNorm( tau );
  
  // Kirchoff pressure. Beware of the sign.
  // Note that this->KSnp1 has been updated by the method
  // KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose )
  double pi = - Trace( this->KSnp1 ) / 3.0;
  
  double oneoverd = 1.0 / this->d( this->pcnp1 );
  double oneovergtaubar = 1.0 / this->g_tau( this->pcnp1 );
  double oneoverm = 1.0 / this->Parameters->flr_m ;
  double oneovergpibar = 1.0 / this->g_pi( pi, this->pcnp1 );
  
  
  //! 3 - Calculate isochoric (blocks 2.1, 2.2) and volumetric (blocks 2.3, 2.4) blocks
  //!     as well as the right hand sides
  
  // 2.1 and 2.2
  ttl::Tensor<4, dim, double> Blockiso = {};
  ttl::Tensor<4, dim, double> RHSiso = {};
  
  FTensors devS = dev( this->Snp1 );
  double normPsid = FrobeniusNorm( devS );
  
  if (
      ( std::fpclassify( taubar ) == FP_ZERO )
      ||
      ( std::fpclassify( normPsid ) == FP_ZERO )
      )
  {
    // if the deviatoric part of S is zero, or if the deviatoric part of KS is zero,
    // the block iso should be zero as well
  }
  else
  {
    double oneovertaubar = 1.0 / taubar;

    // Zwxcoeff1 is the coefficient of Zwx at page 19 of the consistent tangent operator notes,
    // in the block 2.1. It has been studied at page 26 of the same notes
    double Zwxcoeff1 = (
                        oneoverd * oneoverd * this->DdDpc( this->pcnp1 ) -
                        ( 1 - oneoverd ) * oneoverm * oneovergtaubar * this->Dg_tauDpc( this->pcnp1 )
                        )
    * pow( taubar * oneovergtaubar, oneoverm ) ;
    
    // Majorcoeff1 is the coefficient of Major_cdwx at page 19 of the consistent tangent operator notes,
    // in the block 2.1.
    double Majorcoeff1 = ( 1 - oneoverd ) * oneoverm * pow( taubar * oneovergtaubar, oneoverm - 1 ) * oneovergtaubar ;
    
    Blockiso( i,j,w,x ) = ( Zwxcoeff1 * Zwx( w,x ) + Majorcoeff1 * tau( k,l ) * Majorcdwx( k,l,w,x ) * oneovertaubar ) * devS( i,j ) * this->Parameters->flr_gamma0 / normPsid;
    Blockiso( i,j,w,x ) = Blockiso( i,j,w,x ) + this->gammadot_d( taubar, this->pcnp1) * Phiabzy( i,j,z,y ) * Lambdazywx( z,y,w,x );
    
    RHSiso( i,j,w,x ) = Majorcoeff1 * tau( k,l ) * Minorcdkl( k,l,w,x ) * oneovertaubar * devS( i,j ) * this->Parameters->flr_gamma0 / normPsid ;
    RHSiso( i,j,w,x ) = RHSiso( i,j,w,x ) + this->gammadot_d( taubar, this->pcnp1) * Phiabzy( i,j,z,y ) * Vzykl( z,y,w,x );
    
  }
  
  
  // 2.3 and 2.4
  ttl::Tensor<4, dim, double> Blockvol = {};
  ttl::Tensor<4, dim, double> RHSvol = {};
  
  if (  ( pi < this->pi_m( this->pcnp1 ) ) || ( (this->pcnp1) >= this->Parameters->cf_pcinf ) )
  { }
  else
  {
    // Zwxcoeff3 is the coefficient of Zwx at page 19 of the consistent tangent operator notes,
    // in the block 2.3. It has been studied at page 26 of the same notes
    double Zwxcoeff3 = - ( this->Dpi_mDpc( this->pcnp1 ) + ( pi - this->pi_m( this->pcnp1 ) ) * oneovergpibar * this->Dg_piDpc( pi, this->pcnp1 ) )  * oneovergpibar ;
    
    // Commoncoeff3 is defined at page 19 of the consistent tangent operator notes,
    // in the block 2.3.
    double Commoncoeff3 =  pow( ( pi - this->pi_m( this->pcnp1 ) ) * oneovergpibar , oneoverm - 1 ) * oneoverm ;
    
    FTensors Psiv =  - onethird * this->betaC( this->pcnp1 ) * ttl::identity<dim>( i,j );
    
    Blockvol( i,j,w,x ) = this->Parameters->flr_gamma0 * Psiv( i,j ) * Commoncoeff3 * ( Zwxcoeff3 * Zwx( w,x ) - onethird * oneovergpibar * ttl::identity<dim>( k,l ) * Majorcdwx( k,l,w,x )  );
    Blockvol( i,j,w,x ) = Blockvol( i,j,w,x ) + onethird * this->gammadot_v( pi, this->pcnp1) *  this->Parameters->cf_g0 / this->Parameters->cf_pcinf * ttl::identity<dim>( i,j ) * Zwx( w,x ) ;
    
    RHSvol( i,j,w,x ) = this->Parameters->flr_gamma0 * Psiv( i,j ) * Commoncoeff3 * ( - onethird * oneovergpibar * ttl::identity<dim>( k,l ) * Minorcdkl( k,l,w,x ) );
  }
  
  
  //! 4 - Solve the system iteratively and return
  
  ttl::Tensor<4, dim, double> formerdmdf = initdmdf;
  ttl::Tensor<4, dim, double> dmdf = {};
  ttl::Tensor<4, dim, double> dmdfdiff = {};
  
  unsigned MaxIt = 500;
  double ItTol = 1e-8;
  
  for ( unsigned icount = 0; icount< MaxIt; icount++ )
  {
    dmdf( i,j,w,x ) = - dt * ( Blockiso( i,j,k,l ) + Blockvol( i,j,k,l ) ) * formerdmdf( k,l,w,x ) - dt * ( RHSiso( i,j,w,x ) + RHSvol( i,j,w,x ) );
    
    dmdfdiff( i,j,w,x ) = dmdf( i,j,w,x ) - formerdmdf( i,j,w,x );
    
    double normcheck  = FrobeniusNorm( dmdfdiff );
    
    if ( Verbose )
    {
      std::cout << " It = " << icount << ", FrobeniusNorm( dmdf - formerdmdf ) = " << normcheck << "\n";
    }
    
    
    if ( normcheck < ItTol )
      break;
    
    formerdmdf = dmdf;
    
  }
  
  
  return dmdf;
}








template <int dim>
void KMS_IJSS2017_Implicit_BE_Staggered<dim>::ttlinverseExceptionHandling( const ttl::Tensor<4, dim, double>& ttlsysmat )
//! Exception handler for the inversion of the fourth order tensor during the evaluation of DMDF
//!
{
  
  //this->IO->log() << "\n\n  ----------------------------- \n";
  //this->IO->log() << "  ttl::inverse Exception Handling \n";
  //this->IO->log() << "  the code is NOT aborting at \n\n  ";
  //this->IO->log() << "  timestep = " <<  this->TimeIntegrationData->TimeStep() << ", time = " << this->TimeIntegrationData->CurrentTime() ;
  //this->IO->log() << "\n\n  but lapack cannot invert the fourth order tensor during the evaluation of DMDF \n";
  //this->IO->log() << "  within the Consistent Tangent Stiffness Matrix for the  \n";
  //this->IO->log() << "  KMS_IJSS2017_Implicit_BE_Staggered integrator. Here is the tensor in a matrix form\n";
  
//  typedef mtl::mat::dense2D<double> Matrix;
//  // typedef mtl::dense_vector<double> Vector;
//  
//  // Matrix analysis - using mtl
//  Matrix SysMat( dim*dim, dim*dim );
//  
//  unsigned r=0;
//  for (unsigned i=0; i<dim*dim; i++)
//    for (unsigned j=0; j<dim*dim; j++)
//      SysMat[i][j] = ttlsysmat.data[r++];
//  
//  this->IO->log() << "\n  system matrix: \n" << SysMat << "\n";
//  
//  // Eigenvalues analysis - using mtl
//  mtl::mat::eigenvalue_solver<Matrix> E1( SysMat );
//  E1.setMaxIteration(100);
//  E1.setTolerance(1.0e-10);
//  E1.calc();
//  this->IO->log() << " eigenvalues: " << E1.get_eigenvalues() << "\n";
  
  //this->IO->log() << "\n\n  ----------------------------- \n\n\n ";
}



template <int dim>
ttl::Tensor<4, dim, double> KMS_IJSS2017_Implicit_BE_Staggered<dim>::DMDF(
                                                                const double dt,
                                                                const bool Verbose
                                                                )
//! This function evaluates the fourth order tensor DM/DF for the KMS_IJSS2017 model
//! This method expects that the Step has already been updated elsewhere, since it uses quantities at step n+1
//! that MUST have been updated before
//!
//! 1 - Calculate intermediate operators
//! 2 - Calculate some derivatives wrt pc
//! 3 - Calculate blocks 2.1, 2.2, 2.3, 2.4 and right hand sides
//! 4 - Solve the linear system and return
//!
//! According to "Debugging and code optimization" notes, page 5, operators Zwx, Major, Lambda
//! have been hatted, and the amount dHpcdpc factorized.
//! This added step 0. Evaluation of dHpcdpc
//!
{
  
  
  using namespace ttlindexes;
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'y'> y;
  
  double onethird = 1.0 / 3.0;
  
  if ( Verbose )
  {
    std::cout << "\n\n  Debugging KMS_IJSS2017_Implicit_BE_Staggered<dim>::DMDF \n";
    
    
    std::string Fstr="";
    
    PrintInVectorForm( this->Fnp1 , Fstr );
    std::cout << "  Fnp1 = " << Fstr;
    Fstr="";
    
    PrintInVectorForm( this->eFnp1 , Fstr );
    std::cout << "  eFnp1 = " << Fstr;
    Fstr="";
    
    PrintInVectorForm( this->pFnp1 , Fstr );
    std::cout << "  pFnp1 = " << Fstr;
    Fstr="";
    
    PrintInVectorForm( this->Snp1 , Fstr );
    std::cout << "  Snp1 = " << Fstr;
    Fstr="";
    
  }
    
  // 0 - Evaluation of dHpcdpc
  
  // Material parameters
  double a1=this->Parameters->hr_a1, a2=this->Parameters->hr_a2, l1=this->Parameters->hr_Lambda1, l2=this->Parameters->hr_Lambda2;

  // derivative dHpc/dpc
  double dHpcdpc = a1 * l1 * exp( - l1 / this->pcnp1 ) + a2 * l2 * exp( - l2 / this->pcnp1 ) ;
    
    
  // 1 - Calculate intermediate operators
    
  // tensor M = pFn * InvpFnp1
  FTensors InvpFnp1 = ttl::inverse( this->pFnp1 );
  FTensors M = this->pFn( i,k ) * InvpFnp1( k,j );
  
  if ( Verbose )
  {
    std::string Fstr="";
    
    PrintInVectorForm( M , Fstr );
    std::cout << "  M = " << Fstr;
    Fstr="";
  }

  // Zwx has been defined at page 17b of the consistent tangent operator notes
  FTensors HattedZwx = this->UpdateHattedZwx( M );
  
  // Lambdazywx has been defined at page 16 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> HattedLambdazywx;
  ttl::Tensor<4, dim, double> dsdm = this->DSDM( this->Fnp1, this->eFnp1, this->pFnp1, this->pcnp1 );
  HattedLambdazywx( z,y,w,x ) = dHpcdpc * dsdm( z,y,w,x ) + this->DSDpc() ( z,y ) * HattedZwx( w,x );
  
  // Vzykl has been defined at page 16 of the consistent tangent operator notes
  // Note that dSdF is the partial derivative and differs from DSDF
  ttl::Tensor<4, dim, double> Vzykl = this->dSdF( this->Fnp1, this->eFnp1, this->pFnp1, this->pcnp1 );
  
  // Minorcdkl and Majorcdwx have been defined at pages 16-17 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> Minorcdkl = this->Minor_cdkl( this->Fnp1, this->eFnp1, this->pFnp1, this->Snp1, M, Vzykl );
  ttl::Tensor<4, dim, double> HattedMajorcdwx = this->HattedMajor_cdwx( this->Fnp1, this->eFnp1, this->pFnp1, this->Snp1, M, HattedLambdazywx, dHpcdpc );
  
  // Phiabzy has been defined at page 18 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> Phiabzy = this->Phi_abzy( this->Snp1, false );
  
  
  // 2 - Calculate some derivatives wrt pc
  
  // taubar is the Frobenius norm of the deviatoric part of the Kirchoff stress.
  // therefor if it is = 0, the Kirchoff stress is purely volumetric
  // Note that this->KSnp1 has been updated by the method
  // KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose )
  FTensors tau = dev( this->KSnp1 );
  double taubar = FrobeniusNorm( tau );
  
  // Kirchoff pressure. Beware of the sign.
  // Note that this->KSnp1 has been updated by the method
  // KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose )
  double pi = - Trace( this->KSnp1 ) / 3.0;
  
  // these denominators should never be zero.
  double oneoverd = 1.0 / this->d( this->pcnp1 );
  double oneovergtaubar = 1.0 / this->g_tau( this->pcnp1 );
  double oneoverm = 1.0 / this->Parameters->flr_m ;
  double oneovergpibar = 1.0 / this->g_pi( pi, this->pcnp1 );
  
  
  //! 3 - Calculate isochoric (blocks 2.1, 2.2) and volumetric (blocks 2.3, 2.4) blocks
  //!     as well as the right hand sides
  
  // 2.1 and 2.2
  ttl::Tensor<4, dim, double> Blockiso = {};
  ttl::Tensor<4, dim, double> RHSiso = {};
  
  FTensors devS = dev( this->Snp1 );
  double normPsid = FrobeniusNorm( devS );
  
  if (
      ( std::fpclassify( taubar ) == FP_ZERO )
      ||
      ( std::fpclassify( normPsid ) == FP_ZERO )
    )
  {
    // if the deviatoric part of S is zero, or if the deviatoric part of KS is zero,
    // the block iso should be zero as well
  }
  else
  {
    double oneovertaubar = 1.0 / taubar;

    // Zwxcoeff1 is the coefficient of Zwx at page 19 of the consistent tangent operator notes,
    // in the block 2.1. It has been studied at page 26 of the same notes
    double Zwxcoeff1 = (
                        oneoverd * oneoverd * this->DdDpc( this->pcnp1 ) -
                        ( 1 - oneoverd ) * oneoverm * oneovergtaubar * this->Dg_tauDpc( this->pcnp1 )
                        )
    * pow( taubar * oneovergtaubar, oneoverm ) ;
    
    // Majorcoeff1 is the coefficient of Major_cdwx at page 19 of the consistent tangent operator notes,
    // in the block 2.1.
    double Majorcoeff1 = ( 1 - oneoverd ) * oneoverm * pow( taubar * oneovergtaubar, oneoverm - 1 ) * oneovergtaubar ;
    
    Blockiso( i,j,w,x ) = ( Zwxcoeff1 * HattedZwx( w,x ) + Majorcoeff1 * tau( k,l ) * HattedMajorcdwx( k,l,w,x ) * oneovertaubar ) * devS( i,j ) * this->Parameters->flr_gamma0 / normPsid;
    Blockiso( i,j,w,x ) = Blockiso( i,j,w,x ) + this->gammadot_d( taubar, this->pcnp1) * Phiabzy( i,j,z,y ) * HattedLambdazywx( z,y,w,x );
    
    RHSiso( i,j,w,x ) = Majorcoeff1 * tau( k,l ) * Minorcdkl( k,l,w,x ) * oneovertaubar * devS( i,j ) * this->Parameters->flr_gamma0 / normPsid ;
    RHSiso( i,j,w,x ) = RHSiso( i,j,w,x ) + this->gammadot_d( taubar, this->pcnp1) * Phiabzy( i,j,z,y ) * Vzykl( z,y,w,x );
    
  }
  
  
  // 2.3 and 2.4
  ttl::Tensor<4, dim, double> Blockvol = {};
  ttl::Tensor<4, dim, double> RHSvol = {};
  
  double pim = this->pi_m( this->pcnp1 );
    
  if (  ( pi < pim ) || ( (this->pcnp1) >= this->Parameters->cf_pcinf ) )
  { }
  else
  {
    // Zwxcoeff3 is the coefficient of Zwx at page 19 of the consistent tangent operator notes,
    // in the block 2.3. It has been studied at page 26 of the same notes
    double Zwxcoeff3 = - ( this->Dpi_mDpc( this->pcnp1 ) + ( pi - pim ) * oneovergpibar * this->Dg_piDpc( pi, this->pcnp1 ) )  * oneovergpibar ;
    
    // Commoncoeff3 is defined at page 19 of the consistent tangent operator notes,
    // in the block 2.3.
    double Commoncoeff3 =  pow( ( pi - pim ) * oneovergpibar , oneoverm - 1 ) * oneoverm ;
    
    FTensors Psiv =  - onethird * this->betaC( this->pcnp1 ) * ttl::identity<dim>( i,j );
    
    Blockvol( i,j,w,x ) = this->Parameters->flr_gamma0 * Psiv( i,j ) * Commoncoeff3 * ( Zwxcoeff3 * HattedZwx( w,x ) - onethird * oneovergpibar * ttl::identity<dim>( k,l ) * HattedMajorcdwx( k,l,w,x )  );
    Blockvol( i,j,w,x ) = Blockvol( i,j,w,x ) + onethird * this->gammadot_v( pi, this->pcnp1) *  this->Parameters->cf_g0 / this->Parameters->cf_pcinf * ttl::identity<dim>( i,j ) * HattedZwx( w,x ) ;
    
    RHSvol( i,j,w,x ) = this->Parameters->flr_gamma0 * Psiv( i,j ) * Commoncoeff3 * ( - onethird * oneovergpibar * ttl::identity<dim>( k,l ) * Minorcdkl( k,l,w,x ) );
  }
  
  
  //! 4 - Solve the system and return
  
  ttl::Tensor<4, dim, double> sysmat = {};
  ttl::Tensor<4, dim, double> invsysmat = {};
  ttl::Tensor<4, dim, double> rhs = {};
  
  sysmat( i,j,w,x ) = ttl::identity<dim>( i,j,w,x ) * dHpcdpc + dt * ( Blockiso( i,j,w,x ) + Blockvol( i,j,w,x ) );
  
  if ( Verbose )
  {
    // std::cout << "\n  sysmat = " << sysmat( i,j,k,l ) << "\n";
    std::string Fstr="";
    
    PrintInVectorForm( sysmat , Fstr );
    std::cout << "  sysmat = " << Fstr;
    Fstr="";
    std::cout << std::setw(20) << std::setprecision(15) << "  FrobeniusNorm( sysmat ) = " << FrobeniusNorm( sysmat ) << "\n";

  }
  
  try
  {
   invsysmat = ttl::inverse( sysmat );
  }
  catch (int i)
  {
    ttlinverseExceptionHandling( sysmat );
    
    return ttl::zero(  w,x,k,l  );
  }

  rhs( i,j,w,x ) = - dt * dHpcdpc * ( RHSiso( i,j,w,x ) + RHSvol( i,j,w,x ) );
  
  ttl::Tensor<4, dim, double> dmdf;
  dmdf( w,x,k,l ) = invsysmat( w,x,i,j ) * rhs( i,j,k,l ) ;
  
  return dmdf;
}









template <int dim>
void KMS_IJSS2017_Implicit_BE_Staggered<dim>::DMDFandDSDF(
                                                ttl::Tensor<4, dim, double> & dmdf,
                                                ttl::Tensor<4, dim, double> & dsdf,
                                                const double dt,
                                                const bool Verbose
                                                )
  //! This function evaluates the fourth order tensor DM/DF for the KMS_IJSS2017 model
  //! This method expects that the Step has already been updated elsewhere, since it uses quantities at step n+1
  //! that MUST have been updated before
  //!
  //! 1 - Calculate intermediate operators
  //! 2 - Calculate some derivatives wrt pc
  //! 3 - Calculate blocks 2.1, 2.2, 2.3, 2.4 and right hand sides
  //! 4 - Solve the system and calculate dmdf
  //! 5 - Calculate dsdf
{
  
  
  using namespace ttlindexes;
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'y'> y;
  
  double onethird = 1.0 / 3.0;
  
  if ( Verbose )
  {
    std::cout << "\n\n  Debugging KMS_IJSS2017_Implicit_BE_Staggered<dim>::DMDF \n";
    
    
    std::string Fstr="";
    
    PrintInVectorForm( this->Fnp1 , Fstr );
    std::cout << "  Fnp1 = " << Fstr;
    Fstr="";
    
    PrintInVectorForm( this->eFnp1 , Fstr );
    std::cout << "  eFnp1 = " << Fstr;
    Fstr="";
    
    PrintInVectorForm( this->pFnp1 , Fstr );
    std::cout << "  pFnp1 = " << Fstr;
    Fstr="";
    
    PrintInVectorForm( this->Snp1 , Fstr );
    std::cout << "  Snp1 = " << Fstr;
    Fstr="";
    
  }
  
    // 0 - Evaluation of dHpcdpc
  
    // Material parameters
  double a1=this->Parameters->hr_a1, a2=this->Parameters->hr_a2, l1=this->Parameters->hr_Lambda1, l2=this->Parameters->hr_Lambda2;
  
    // derivative dHpc/dpc
  double dHpcdpc = a1 * l1 * exp( - l1 / this->pcnp1 ) + a2 * l2 * exp( - l2 / this->pcnp1 ) ;
  
  
    // 1 - Calculate intermediate operators
  
    // tensor M = pFn * InvpFnp1
  FTensors InvpFnp1 = ttl::inverse( this->pFnp1 );
  FTensors M = this->pFn( i,k ) * InvpFnp1( k,j );
  
  if ( Verbose )
  {
    std::string Fstr="";
    
    PrintInVectorForm( M , Fstr );
    std::cout << "  M = " << Fstr;
    Fstr="";
  }
  
    // Zwx has been defined at page 17b of the consistent tangent operator notes
  FTensors HattedZwx = this->UpdateHattedZwx( M );
  
    // Lambdazywx has been defined at page 16 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> HattedLambdazywx;
  ttl::Tensor<4, dim, double> dsdm = this->DSDM( this->Fnp1, this->eFnp1, this->pFnp1, this->pcnp1 );
  HattedLambdazywx( z,y,w,x ) = dHpcdpc * dsdm( z,y,w,x ) + this->DSDpc() ( z,y ) * HattedZwx( w,x );
  
    // Vzykl has been defined at page 16 of the consistent tangent operator notes
    // Note that dSdF is the partial derivative and differs from DSDF
    // Vzykl can be seen as the "elastic" change in the second Piola-Kirchoff stress, in the sense
    // that if an increment in F is purely elastic (DF->DeF) it causes neither changes in M nor in pc
    // and the stress increment coincides with Vzykl.
  ttl::Tensor<4, dim, double> Vzykl = this->dSdF( this->Fnp1, this->eFnp1, this->pFnp1, this->pcnp1 );
  
    // Minorcdkl and Majorcdwx have been defined at pages 16-17 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> Minorcdkl = this->Minor_cdkl( this->Fnp1, this->eFnp1, this->pFnp1, this->Snp1, M, Vzykl );
  ttl::Tensor<4, dim, double> HattedMajorcdwx = this->HattedMajor_cdwx( this->Fnp1, this->eFnp1, this->pFnp1, this->Snp1, M, HattedLambdazywx, dHpcdpc );
  
    // Phiabzy has been defined at page 18 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> Phiabzy = this->Phi_abzy( this->Snp1, false );
  
  
    // 2 - Calculate some derivatives wrt pc
  
    // taubar is the Frobenius norm of the deviatoric part of the Kirchoff stress.
    // therefor if it is = 0, the Kirchoff stress is purely volumetric
    // Note that this->KSnp1 has been updated by the method
    // KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose )
  FTensors tau = dev( this->KSnp1 );
  double taubar = FrobeniusNorm( tau );
  
    // Kirchoff pressure. Beware of the sign.
    // Note that this->KSnp1 has been updated by the method
    // KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose )
  double pi = - Trace( this->KSnp1 ) / 3.0;
  
    // these denominators should never be zero.
  double oneoverd = 1.0 / this->d( this->pcnp1 );
  double oneovergtaubar = 1.0 / this->g_tau( this->pcnp1 );
  double oneoverm = 1.0 / this->Parameters->flr_m ;
  double oneovergpibar = 1.0 / this->g_pi( pi, this->pcnp1 );
  
  
    //! 3 - Calculate isochoric (blocks 2.1, 2.2) and volumetric (blocks 2.3, 2.4) blocks
    //!     as well as the right hand sides
  
    // 2.1 and 2.2
  ttl::Tensor<4, dim, double> Blockiso = {};
  ttl::Tensor<4, dim, double> RHSiso = {};
  
  FTensors devS = dev( this->Snp1 );
  double normPsid = FrobeniusNorm( devS );
  
  if (
      ( std::fpclassify( taubar ) == FP_ZERO )
      ||
      ( std::fpclassify( normPsid ) == FP_ZERO )
      )
  {
      // if the deviatoric part of S is zero, or if the deviatoric part of KS is zero,
      // the block iso should be zero as well
  }
  else
  {
    double oneovertaubar = 1.0 / taubar;
    
      // Zwxcoeff1 is the coefficient of Zwx at page 19 of the consistent tangent operator notes,
      // in the block 2.1. It has been studied at page 26 of the same notes
    double Zwxcoeff1 = (
                        oneoverd * oneoverd * this->DdDpc( this->pcnp1 ) -
                        ( 1 - oneoverd ) * oneoverm * oneovergtaubar * this->Dg_tauDpc( this->pcnp1 )
                        )
    * pow( taubar * oneovergtaubar, oneoverm ) ;
    
      // Majorcoeff1 is the coefficient of Major_cdwx at page 19 of the consistent tangent operator notes,
      // in the block 2.1.
    double Majorcoeff1 = ( 1 - oneoverd ) * oneoverm * pow( taubar * oneovergtaubar, oneoverm - 1 ) * oneovergtaubar ;
    
    Blockiso( i,j,w,x ) = ( Zwxcoeff1 * HattedZwx( w,x ) + Majorcoeff1 * tau( k,l ) * HattedMajorcdwx( k,l,w,x ) * oneovertaubar ) * devS( i,j ) * this->Parameters->flr_gamma0 / normPsid;
    Blockiso( i,j,w,x ) = Blockiso( i,j,w,x ) + this->gammadot_d( taubar, this->pcnp1) * Phiabzy( i,j,z,y ) * HattedLambdazywx( z,y,w,x );
    
    RHSiso( i,j,w,x ) = Majorcoeff1 * tau( k,l ) * Minorcdkl( k,l,w,x ) * oneovertaubar * devS( i,j ) * this->Parameters->flr_gamma0 / normPsid ;
    RHSiso( i,j,w,x ) = RHSiso( i,j,w,x ) + this->gammadot_d( taubar, this->pcnp1) * Phiabzy( i,j,z,y ) * Vzykl( z,y,w,x );
    
  }
  
  
    // 2.3 and 2.4
  ttl::Tensor<4, dim, double> Blockvol = {};
  ttl::Tensor<4, dim, double> RHSvol = {};
  
  double pim = this->pi_m( this->pcnp1 );
  
  if (  ( pi < pim ) || ( (this->pcnp1) >= this->Parameters->cf_pcinf ) )
  { }
  else
  {
      // Zwxcoeff3 is the coefficient of Zwx at page 19 of the consistent tangent operator notes,
      // in the block 2.3. It has been studied at page 26 of the same notes
    double Zwxcoeff3 = - ( this->Dpi_mDpc( this->pcnp1 ) + ( pi - pim ) * oneovergpibar * this->Dg_piDpc( pi, this->pcnp1 ) )  * oneovergpibar ;
    
      // Commoncoeff3 is defined at page 19 of the consistent tangent operator notes,
      // in the block 2.3.
    double Commoncoeff3 =  pow( ( pi - pim ) * oneovergpibar , oneoverm - 1 ) * oneoverm ;
    
    FTensors Psiv =  - onethird * this->betaC( this->pcnp1 ) * ttl::identity<dim>( i,j );
    
    Blockvol( i,j,w,x ) = this->Parameters->flr_gamma0 * Psiv( i,j ) * Commoncoeff3 * ( Zwxcoeff3 * HattedZwx( w,x ) - onethird * oneovergpibar * ttl::identity<dim>( k,l ) * HattedMajorcdwx( k,l,w,x )  );
    Blockvol( i,j,w,x ) = Blockvol( i,j,w,x ) + onethird * this->gammadot_v( pi, this->pcnp1) *  this->Parameters->cf_g0 / this->Parameters->cf_pcinf * ttl::identity<dim>( i,j ) * HattedZwx( w,x ) ;
    
    RHSvol( i,j,w,x ) = this->Parameters->flr_gamma0 * Psiv( i,j ) * Commoncoeff3 * ( - onethird * oneovergpibar * ttl::identity<dim>( k,l ) * Minorcdkl( k,l,w,x ) );
  }
  
  
    //! 4 - Solve the system and return
  
  ttl::Tensor<4, dim, double> sysmat = {};
  ttl::Tensor<4, dim, double> invsysmat = {};
  ttl::Tensor<4, dim, double> rhs = {};
  
  sysmat( i,j,w,x ) = ttl::identity<dim>( i,j,w,x ) * dHpcdpc + dt * ( Blockiso( i,j,w,x ) + Blockvol( i,j,w,x ) );
  
  if ( Verbose )
  {
      // std::cout << "\n  sysmat = " << sysmat( i,j,k,l ) << "\n";
    std::string Fstr="";
    
    PrintInVectorForm( sysmat , Fstr );
    std::cout << "  sysmat = " << Fstr;
    Fstr="";
    std::cout << std::setw(20) << std::setprecision(15) << "  FrobeniusNorm( sysmat ) = " << FrobeniusNorm( sysmat ) << "\n";
    
  }
  
  try
  {
    invsysmat = ttl::inverse( sysmat );
  }
  catch (int i)
  {
    ttlinverseExceptionHandling( sysmat );
    
    dmdf( w,x,k,l ) = ttl::zero(  w,x,k,l  );
    dsdf( z,y,k,l ) = Vzykl( z,y,k,l );
    return;
  }
  
  rhs( i,j,w,x ) = - dt * dHpcdpc * ( RHSiso( i,j,w,x ) + RHSvol( i,j,w,x ) );
  
  dmdf( w,x,k,l ) = invsysmat( w,x,i,j ) * rhs( i,j,k,l ) ;

  
    //! 5 - Calculate dsdf
  double FrNormdmdf = FrobeniusNorm( dmdf );
  double oneoverdHpcdpc = 1.0 / dHpcdpc;
  
  if ( std::fpclassify( FrNormdmdf ) == FP_ZERO )
  {
    dsdf( z,y,k,l ) = Vzykl( z,y,k,l );
    return;
  }
  else if ( ( FrNormdmdf < 1e-16 ) || ( dHpcdpc <  1e-16 ) )
  {
    
    //this->IO->log() << "\n\n  ----------------------------- \n";
    //this->IO->log() << "  dsdf Exception Handling \n";
    //this->IO->log() << "  the code is NOT aborting at \n\n  ";
    //this->IO->log() << "  timestep = " <<  this->TimeIntegrationData->TimeStep() << ", time = " << this->TimeIntegrationData->CurrentTime() ;
    //this->IO->log() << "\n\n  but outcomes might be severely wrong because ";
    //this->IO->log() << "  ( FrNormdmdf < 1e-16 ) || ( dHpcdpc <  1e-16 ) during the evaluation of dsdf \n";
    //this->IO->log() << "  within the Consistent Tangent Stiffness Matrix for the  \n";
    //this->IO->log() << "  KMS_IJSS2017_Implicit_BE_Staggered integrator.\n";
    
    if ( ( FrNormdmdf < 1e-16 ) && ( this->pcnp1 < this->Parameters->d_pcb ) )
    {
      dsdf( z,y,k,l ) = Vzykl( z,y,k,l );

      //this->IO->log() << "  Since pc < pcb and FrNormdmdf < 1e-16, it has been arbitrarily taken   \n";
      //this->IO->log() << "  that the change in the second Piola-Kirchoff stress is elastic in nature, i.e.\n";
      //this->IO->log() << "  dsdf( z,y,k,l ) = Vzykl( z,y,k,l )   \n";

    }
    else
    {
      dsdf( z,y,k,l ) = oneoverdHpcdpc * HattedLambdazywx( z,y,w,x ) * dmdf( w,x,k,l ) + Vzykl( z,y,k,l );
    }
    //this->IO->log() << "\n\n  ----------------------------- \n\n\n ";
    
  }
  else
  {
    dsdf( z,y,k,l ) = oneoverdHpcdpc * HattedLambdazywx( z,y,w,x ) * dmdf( w,x,k,l ) + Vzykl( z,y,k,l );
  }
  
}


template <int dim>
int KMS_IJSS2017_Implicit_BE_Staggered<dim>::compute_dMdF(double *dMdF_in,
                                                          const double dt)
  //! This function evaluates the fourth order tensor dM/Du for the KMS_IJSS2017 model
  //! This method expects that the Step has already been updated elsewhere, since it uses quantities at step n+1
  //! that MUST have been updated before
  //!
  //! 1 - Calculate intermediate operators
  //! 2 - Calculate some derivatives wrt pc
  //! 3 - Calculate blocks 2.1, 2.2, 2.3, 2.4 and right hand sides
  //! 4 - Solve the system and calculate dMdF
{
  int err = 0;
  using namespace ttlindexes;
  static constexpr ttl::Index<'w'> w;
  static constexpr ttl::Index<'x'> x;
  static constexpr ttl::Index<'z'> z;
  static constexpr ttl::Index<'y'> y;
 
  ttl::Tensor<4, dim, double*> dMdF(dMdF_in);
       
  double onethird = 1.0/3.0;
  
    // 0 - Evaluation of dHpcdpc
  
    // Material parameters
  double a1=this->Parameters->hr_a1;
  double a2=this->Parameters->hr_a2;
  double l1=this->Parameters->hr_Lambda1;
  double l2=this->Parameters->hr_Lambda2;
  
  // derivative dHpc/dpc
  double dHpcdpc = a1*l1*exp(-l1/this->pcnp1) + a2*l2*exp(-l2/this->pcnp1) ;
  
  
  // 1 - Calculate intermediate operators  
  // tensor M = pFn * pFnp1I
  FTensors pFnp1I = ttl::inverse(this->pFnp1);
    
  FTensors M = this->pFn(i,k)*pFnp1I(k,j);
  
  // Zwx has been defined at page 17b of the consistent tangent operator notes
  FTensors HattedZwx = this->UpdateHattedZwx(M);
  
  // Lambdazywx has been defined at page 16 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> HattedLambdazywx;
  ttl::Tensor<4, dim, double> dsdm = this->DSDM( this->Fnp1, this->eFnp1, this->pFnp1, this->pcnp1 );
  HattedLambdazywx( z,y,w,x ) = dHpcdpc * dsdm( z,y,w,x ) + this->DSDpc() ( z,y ) * HattedZwx( w,x );
  
  // Vzykl has been defined at page 16 of the consistent tangent operator notes
  // Note that dSdF is the partial derivative and differs from DSDF
  // Vzykl can be seen as the "elastic" change in the second Piola-Kirchoff stress, in the sense
  // that if an increment in F is purely elastic (DF->DeF) it causes neither changes in M nor in pc
  // and the stress increment coincides with Vzykl.
  ttl::Tensor<4, dim, double> Vzykl = this->dSdF( this->Fnp1, this->eFnp1, this->pFnp1, this->pcnp1 );
  
  // Minorcdkl and Majorcdwx have been defined at pages 16-17 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> Minorcdkl = this->Minor_cdkl( this->Fnp1, this->eFnp1, this->pFnp1, this->Snp1, M, Vzykl );
  ttl::Tensor<4, dim, double> HattedMajorcdwx = this->HattedMajor_cdwx( this->Fnp1, this->eFnp1, this->pFnp1, this->Snp1, M, HattedLambdazywx, dHpcdpc );
  
  // Phiabzy has been defined at page 18 of the consistent tangent operator notes
  ttl::Tensor<4, dim, double> Phiabzy = this->Phi_abzy( this->Snp1, false );
  
  
  // 2 - Calculate some derivatives wrt pc
  
  // taubar is the Frobenius norm of the deviatoric part of the Kirchoff stress.
  // therefor if it is = 0, the Kirchoff stress is purely volumetric
  // Note that this->KSnp1 has been updated by the method
  // KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose )
  FTensors tau = dev( this->KSnp1 );
  double taubar = FrobeniusNorm( tau );
  
  // Kirchoff pressure. Beware of the sign.
  // Note that this->KSnp1 has been updated by the method
  // KMS_IJSS2017_Implicit_BE_Staggered<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose )
  double pi = - Trace( this->KSnp1 ) / 3.0;
  
    // these denominators should never be zero.
  double oneoverd = 1.0 / this->d( this->pcnp1 );
  double oneovergtaubar = 1.0 / this->g_tau( this->pcnp1 );
  double oneoverm = 1.0 / this->Parameters->flr_m ;
  double oneovergpibar = 1.0 / this->g_pi( pi, this->pcnp1 );
  
  
  //! 3 - Calculate isochoric (blocks 2.1, 2.2) and volumetric (blocks 2.3, 2.4) blocks
  //!     as well as the right hand sides
  
  // 2.1 and 2.2
  ttl::Tensor<4, dim, double> Blockiso = {};
  ttl::Tensor<4, dim, double> RHSiso = {};
  
  FTensors devS = dev( this->Snp1 );
  double normPsid = FrobeniusNorm( devS );
  
  if (
      ( std::fpclassify( taubar ) == FP_ZERO )
      ||
      ( std::fpclassify( normPsid ) == FP_ZERO )
      )
  {
      // if the deviatoric part of S is zero, or if the deviatoric part of KS is zero,
      // the block iso should be zero as well
  }
  else
  {
    double oneovertaubar = 1.0 / taubar;
    
      // Zwxcoeff1 is the coefficient of Zwx at page 19 of the consistent tangent operator notes,
      // in the block 2.1. It has been studied at page 26 of the same notes
    double Zwxcoeff1 = (
                        oneoverd * oneoverd * this->DdDpc( this->pcnp1 ) -
                        ( 1 - oneoverd ) * oneoverm * oneovergtaubar * this->Dg_tauDpc( this->pcnp1 )
                        )
    * pow( taubar * oneovergtaubar, oneoverm ) ;
    
      // Majorcoeff1 is the coefficient of Major_cdwx at page 19 of the consistent tangent operator notes,
      // in the block 2.1.
    double Majorcoeff1 = ( 1 - oneoverd ) * oneoverm * pow( taubar * oneovergtaubar, oneoverm - 1 ) * oneovergtaubar ;
    
    Blockiso( i,j,w,x ) = ( Zwxcoeff1 * HattedZwx( w,x ) + Majorcoeff1 * tau( k,l ) * HattedMajorcdwx( k,l,w,x ) * oneovertaubar ) * devS( i,j ) * this->Parameters->flr_gamma0 / normPsid;
    Blockiso( i,j,w,x ) = Blockiso( i,j,w,x ) + this->gammadot_d( taubar, this->pcnp1) * Phiabzy( i,j,z,y ) * HattedLambdazywx( z,y,w,x );
    
    RHSiso( i,j,w,x ) = Majorcoeff1 * tau( k,l ) * Minorcdkl( k,l,w,x ) * oneovertaubar * devS( i,j ) * this->Parameters->flr_gamma0 / normPsid ;
    RHSiso( i,j,w,x ) = RHSiso( i,j,w,x ) + this->gammadot_d( taubar, this->pcnp1) * Phiabzy( i,j,z,y ) * Vzykl( z,y,w,x );
    
  }
  
  
    // 2.3 and 2.4
  ttl::Tensor<4, dim, double> Blockvol = {};
  ttl::Tensor<4, dim, double> RHSvol = {};
  
  double pim = this->pi_m( this->pcnp1 );
  
  if (  ( pi < pim ) || ( (this->pcnp1) >= this->Parameters->cf_pcinf ) )
  { }
  else
  {
      // Zwxcoeff3 is the coefficient of Zwx at page 19 of the consistent tangent operator notes,
      // in the block 2.3. It has been studied at page 26 of the same notes
    double Zwxcoeff3 = - ( this->Dpi_mDpc( this->pcnp1 ) + ( pi - pim ) * oneovergpibar * this->Dg_piDpc( pi, this->pcnp1 ) )  * oneovergpibar ;
    
      // Commoncoeff3 is defined at page 19 of the consistent tangent operator notes,
      // in the block 2.3.
    double Commoncoeff3 =  pow( ( pi - pim ) * oneovergpibar , oneoverm - 1 ) * oneoverm ;
    
    FTensors Psiv =  - onethird * this->betaC( this->pcnp1 ) * ttl::identity<dim>( i,j );
    
    Blockvol( i,j,w,x ) = this->Parameters->flr_gamma0 * Psiv( i,j ) * Commoncoeff3 * ( Zwxcoeff3 * HattedZwx( w,x ) - onethird * oneovergpibar * ttl::identity<dim>( k,l ) * HattedMajorcdwx( k,l,w,x )  );
    Blockvol( i,j,w,x ) = Blockvol( i,j,w,x ) + onethird * this->gammadot_v( pi, this->pcnp1) *  this->Parameters->cf_g0 / this->Parameters->cf_pcinf * ttl::identity<dim>( i,j ) * HattedZwx( w,x ) ;
    
    RHSvol( i,j,w,x ) = this->Parameters->flr_gamma0 * Psiv( i,j ) * Commoncoeff3 * ( - onethird * oneovergpibar * ttl::identity<dim>( k,l ) * Minorcdkl( k,l,w,x ) );
  }
  
  
    //! 4 - Solve the system and return
  
  ttl::Tensor<4, dim, double> sysmat = {};
  ttl::Tensor<4, dim, double> invsysmat = {};
  ttl::Tensor<4, dim, double> rhs = {};
  
  sysmat( i,j,w,x ) = ttl::identity<dim>( i,j,w,x ) * dHpcdpc + dt * ( Blockiso( i,j,w,x ) + Blockvol( i,j,w,x ) );
  
  try
  {
    invsysmat = ttl::inverse( sysmat );
  }
  catch (int i)
  {
    //ttlinverseExceptionHandling( sysmat );
    
    dMdF( w,x,k,l ) = ttl::zero(  w,x,k,l  );
    err++;
    return err;
  }
  
  rhs( i,j,w,x ) = - dt * dHpcdpc * ( RHSiso( i,j,w,x ) + RHSvol( i,j,w,x ) );  
  dMdF( w,x,k,l ) = invsysmat( w,x,i,j ) * rhs( i,j,k,l ) ;
  
  return err;
}

template <int dim>
int KMS_IJSS2017_Implicit_BE_Staggered<dim>::update_elasticity(double *eF_in,
                                                               double pc,
                                                               double *eS_in,
                                                               double *L_in,
                                                               const int compute_elasticity)
{  
  int err = 0;

  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;

  // compute PK2 stress
  ttl::Tensor<2, dim, double*> eF(eF_in);
  ttl::Tensor<2, dim, double*> eS(eS_in);

  // Je at step n+1
  double eJ = det(eF);
  
  // Left Cauchy-Green tensor and its inverse
  ttl::Tensor<2, dim, double> C, CI;
  C(i,j) = eF(k,i)*eF(k,j);
  
  try{
    CI = ttl::inverse(C);
  }
  catch(int inv_err)
  {
    err++;
    return err;
  }
  
  update_elasticity_dev(eF.data, pc, eS_in, L_in, compute_elasticity);
  double du  = compute_dudj(eJ, pc);
  double ddu = compute_d2udj2(eJ, pc);
    
  eS(i,j) += du*eJ*CI(i,j) ;
  
  if(compute_elasticity == 0)
    return err;

  ttl::Tensor<4, dim, double*> L(L_in);  
  ttl::Tensor<4, dim, double> CIoxCI, CICI;
    
  CIoxCI(i,j,k,l) = CI(i,j)*CI(k,l);
  CICI(i,j,k,l)   = CI(i,k)*CI(l,j); 
  
  L(i,j,k,l) += (eJ*du + eJ*eJ*ddu)*CIoxCI(i,j,k,l)
                - 2.0*eJ*du*CICI(i,j,k,l);         
  return err;
}                                                               

template <int dim>
int KMS_IJSS2017_Implicit_BE_Staggered<dim>::update_elasticity_dev(double *eF_in,
                                                                   double pc,
                                                                   double *eS_in,
                                                                   double *L_in,
                                                                   const int compute_elasticity)
{  
  int err = 0;

  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'k'> k;
  static constexpr ttl::Index<'l'> l;
  static constexpr ttl::Index<'m'> m;
  static constexpr ttl::Index<'n'> n;
  static constexpr ttl::Index<'o'> o;

  // compute PK2 stress
  ttl::Tensor<2, dim, double*> eF(eF_in);
  this->eFnp1 = eF;
  this->pcnp1 = pc;
      
  ttl::Tensor<2, dim, double*> eS(eS_in);

  FTensors I = ttl::identity(i,j);

  // Je at step n+1
  double eJ = det(eF);
  
  // Left Cauchy-Green tensor and its inverse
  FTensors Ce;
  Ce(i,j) = eF(k,i)*eF(k,j);
  
  FTensors InvCe;
  try{
    InvCe = ttl::inverse(Ce);
  }
  catch(int inv_err)
  {
    err++;
    return err;
  }
  
  double mu    = this->shearmodulus(pc);
  double trCe  = Trace(Ce);
  
  // Isochoric contribution
  eS(i,j) = mu*pow(eJ, -2.0/3.0)*(I(i,j) - trCe/3.0*InvCe(i,j));
  
  if(compute_elasticity == 0)
    return err;
    
  ttl::Tensor<4, dim, double*> dSdC(L_in); 
               
  double onethird = 1.0/3.0;
  double powJe = pow(eJ, -2.0*onethird);
  
  dSdC(i,j,k,l) = 2.0*mu*(-onethird*powJe*InvCe(k,l))*(I(i,j) - trCe*onethird*InvCe(i,j)) +
                  2.0*mu*powJe*onethird*(-I(k,l)*InvCe(i,j) + trCe*InvCe(i,k)*InvCe(l,j));
            
  return err;
} 

template <int dim>
double KMS_IJSS2017_Implicit_BE_Staggered<dim>::compute_dudj(double eJ,
                                                             double pc)
{  
  KMS_IJSS2017_Parameters *P = this->Parameters;
  double kappa = this->bulkmodulus(pc);
  double c_of_pc = this->c(pc); 

  double exponent = - pow((1.0-this->c(pc)/P->c_inf), P->pl_n)/P->K_kappa;  
  double dudj = (0.5*kappa*(eJ*log(eJ)+eJ - 1.0 ) + (c_of_pc-(P->K_p0 + c_of_pc)*pow(eJ,exponent)))/eJ;
  
  return dudj;
}

template <int dim>
double KMS_IJSS2017_Implicit_BE_Staggered<dim>::compute_d2udj2(double eJ,
                                                               double pc)
{  
  KMS_IJSS2017_Parameters *P = this->Parameters;
  double kappa = this->bulkmodulus(pc);
  double c_of_pc = this->c(pc); 
               
  double alpha = pow(1.0 - c_of_pc/P->c_inf, P->pl_n)/P->K_kappa;  
  double Upp = 0.5*kappa/eJ*(1.0+1.0/eJ) - 1.0/eJ/eJ*(c_of_pc-(P->K_p0 + c_of_pc)*(1.0+alpha)*pow(eJ, -alpha));
  
  return Upp;
}
#endif /* KMS_IJSS2017_implicit_staggered_h */
