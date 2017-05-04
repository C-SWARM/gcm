//
//  KMS-IJSS2017-integrator-templates.hpp
//  ttl-learning
//
//  Created by alberto salvadori on 1/18/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

/*!
 
 \file KMS-IJSS2017-integrator-templates.h
 \brief Template code file for the class base of a generic integration algorithm for the KMS-IJSS2017 model
 
 This file contains template classes and methods that pertain to every integrator of the poroviscoplastic model desccribed in the paper
 "A. Krairi, K. Matouš, A. Salvadori, A poro-viscoplastic constitutive model for granular materials at finite strain, submitted to IJSS (2017)".
 
 */


#ifndef KMS_IJSS2017_integrator_templates_h
#define KMS_IJSS2017_integrator_templates_h


#include <iostream>
#include <iomanip>
#include <cstring>





// class KMS_IJSS2017_Integrator
// *****************************


// Constructors
// ------------

template <int dim>
KMS_IJSS2017_Integrator<dim>::KMS_IJSS2017_Integrator()
//!< Default Constructor
{
  IO= NULL;
  pcEQpc_inf = false;
}


template <int dim>
KMS_IJSS2017_Integrator<dim>::KMS_IJSS2017_Integrator( KMS_IJSS2017_IO *InpOut, TimeIntegrationDataManager* TIDM  )
//! Constructor with material parameters
{
  
  // Assignments
  IO = InpOut;
  this->TimeIntegrationData = TIDM;

  pcEQpc_inf = false;
  
  // Definitions
  if ( IO != NULL )
    Print = true;
  else
    Print = false;
  
}


// Public Methods
// --------------



template <int dim>
void KMS_IJSS2017_Integrator<dim>::StepClose( )
//void KMS_IJSS2017_Integrator<dim>::StepClose( )
//! Assigns data at step n from data at step n+1
{
  
  // Internal variable
  pcn = pcnp1;

  // "Global" variables
  FeFp_Integrator<dim>::StepClose( );

}


template <int dim>
void KMS_IJSS2017_Integrator<dim>::PrintData( size_t step, double time, bool PrintOnScreen )
//! Prints data of each step on the proper files
{
  
  std::string str="";
  
  std::ostream& Foutstream = PrintOnScreen ? std::cout : IO->F();
  std::ostream& pFoutstream = PrintOnScreen ? std::cout : IO->pF();
  std::ostream& SPKFoutstream = PrintOnScreen ? std::cout : IO->SPK();
  std::ostream& KSoutstream = PrintOnScreen ? std::cout : IO->KS();
  std::ostream& sigmaoutstream = PrintOnScreen ? std::cout : IO->sigma();
  std::ostream& pcFoutstream = PrintOnScreen ? std::cout : IO->pc();

  
  if ( Print && step == 0)
  {
    // Output for F
    Foutstream << "step, time, J, F(11), F(12), F(13), F(21), F(22), F(23), F(31), F(32), F(33) \n";
    str="";
    PrintInVectorForm( this->Fn, str );
    Foutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << std::setw(20) << std::setprecision(15) << det( this->Fn ) << ", " << str;
    Foutstream.flush();
    
    // Output for Fp
    pFoutstream << "step, time, Jp, Fp(11), Fp(12), Fp(13), Fp(21), Fp(22), Fp(23), Fp(31), Fp(32), Fp(33) \n";
    str="";
    PrintInVectorForm( this->pFn, str );
    pFoutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << std::setw(20) << std::setprecision(15) << det( this->pFn ) << ", " << str;
    pFoutstream.flush();

    
    // Output for Fe if PrintOnScreen
    if ( PrintOnScreen )
    {
      std::cout << "step, time, Je, Fe(11), Fe(12), Fe(13), Fe(21), Fe(22), Fe(23), Fe(31), Fe(32), Fe(33) \n";
      str="";
      PrintInVectorForm( this->eFn, str );
      std::cout << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << std::setw(20) << std::setprecision(15) << det( this->eFn ) << ", " << str;
      std::cout.flush();
    }

    // Output for S
    SPKFoutstream << "step, time, SecondPiolaKirchoff(11), SecondPiolaKirchoff(12), SecondPiolaKirchoff(13), SecondPiolaKirchoff(21), SecondPiolaKirchoff(22), SecondPiolaKirchoff(23), SecondPiolaKirchoff(31), SecondPiolaKirchoff(32), SecondPiolaKirchoff(33) \n";
    str="";
    PrintInVectorForm( this->Sn, str );
    SPKFoutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << str;
    SPKFoutstream.flush();
    
    // Output for KS
    KSoutstream << "step, time, KirchoffStress(11), KirchoffStress(12), KirchoffStress(13), KirchoffStress(21), KirchoffStress(22), KirchoffStress(23), KirchoffStress(31), KirchoffStress(32), KirchoffStress(33) \n";
    str="";
    PrintInVectorForm( this->KSn, str );
    KSoutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << str;
    KSoutstream.flush();

    // Output for sigma
    sigmaoutstream << "step, time, sigma(11), sigma(12), sigma(13), sigma(21), sigma(22), sigma(23), sigma(31), sigma(32), sigma(33) \n";
    str="";
    PrintInVectorForm( this->sigman, str );
    sigmaoutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << str;
    sigmaoutstream.flush();
    
    // Output for pc
    pcFoutstream << "step, time, pc\n";
    pcFoutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << std::setw(20) << std::setprecision(15) << this->pcn << "\n";
    pcFoutstream.flush();

    
  }
  
  
  else if ( Print && IO->StepToPrint(step) )
  {
    // Output for F
    str="";
    PrintInVectorForm( this->Fn, str );
    Foutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << std::setw(20) << std::setprecision(15) << det( this->Fn )  << ", " << str;
    Foutstream.flush();
    
    // Output for Fp
    str="";
    PrintInVectorForm( this->pFn, str );
    pFoutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << std::setw(20) << std::setprecision(15) << det( this->pFn )  << ", " << str;
    pFoutstream.flush();

    // Output for S
    str="";
    PrintInVectorForm( this->Sn, str );
    SPKFoutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << str;
    SPKFoutstream.flush();
    
    // Output for KS
    str="";
    PrintInVectorForm( this->KSn, str );
    KSoutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << str;
    KSoutstream.flush();
    
    // Output for sigma
    str="";
    PrintInVectorForm( this->sigman, str );
    sigmaoutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << str;
    sigmaoutstream.flush();

    // Output for pc
    pcFoutstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << std::setw(20) << std::setprecision(15) << this->pcn << "\n";
    pcFoutstream.flush();

  }
  
}



// class KMS_IJSS2017_Integration_Algorithms
// *****************************************


// Constructors
// ------------

template <int dim>
KMS_IJSS2017_Integration_Algorithms<dim>::KMS_IJSS2017_Integration_Algorithms(KMS_IJSS2017_Parameters* P, KMS_IJSS2017_IO* InpOut, TimeIntegrationDataManager* TIDM )
 : KMS_IJSS2017<dim> ( P ), KMS_IJSS2017_Integrator<dim> ( InpOut, TIDM )
{  
}
//! Constructor with material parameters and IO


template <int dim>
void KMS_IJSS2017_Integration_Algorithms<dim>::PrintStep( size_t step, double time, bool PrintOnScreen )
//! Prints data and material properties of each step on the proper files
{
  this->PrintData( step, time, PrintOnScreen );
  
  std::ostream& outstream = PrintOnScreen ? std::cout : this->IO->mat();
  
  if ( this->Print && step == 0)
  {
    // Output for material properties
    outstream << "step, time, pc, d, c, mu, k, \n";
    outstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << std::setw(20) << std::setprecision(15) << this->pcn
                            << ", " << std::setw(20) << std::setprecision(15) << this->d(this->pcn) << ", " << std::setw(20) << std::setprecision(15) << this->c(this->pcn)
    << ", " << std::setw(20) << std::setprecision(15) << this->shearmodulus(this->pcn) << ", " << std::setw(20) << std::setprecision(15) << this->bulkmodulus(this->pcn) << "\n";
    outstream.flush();
  }
  else if ( this->Print && this->IO->StepToPrint(step) )
  {
    // Output for pc
    outstream << step << ", " << std::setw(20) << std::setprecision(15) << time << ", " << std::setw(20) << std::setprecision(15) << this->pcn
                            << ", " << std::setw(20) << std::setprecision(15) << this->d(this->pcn) << ", " << std::setw(20) << std::setprecision(15) << this->c(this->pcn)
    << ", " << std::setw(20) << std::setprecision(15) << this->shearmodulus(this->pcn) << ", " << std::setw(20) << std::setprecision(15) << this->bulkmodulus(this->pcn) << "\n";
    outstream.flush();

  }

}



// Generic methods
// ---------------


template <int dim>
void KMS_IJSS2017_Integration_Algorithms<dim>::Initialize( const double pc0, const FTensors& F0, const FTensors& pF0, const FTensors& S0  )
//! Initialize the integrator with pcn, Fn, and Sn at t=0
//! The quantities at n+1 are also initialized pcn, Fn, and Sn at t=0
//! so that the consistent tangent stiffness at the very initial time can be estimated
//! as usual
{

  using namespace ttlindexes;
  
  // Zeroing
  this->Dpn = ttl::identity(i,j);
  
  // Assignment of pcn
  this->pcn = pc0;
  this->pcnp1 = pc0;
  
  // Assignment of F tensors
  this->Fn = F0;
  this->Fnp1 = F0;
  this->pFn = pF0;
  this->pFnp1 = pF0;
  
  // Check the compatibility for the plastic tensor
  bool JpInitWarning = false;
  double Jperr = fabs( log( det( pF0 ) ) - this->HardeningLaw( pc0 ) );
  if ( Jperr / fabs( log( det( pF0 ) ) ) > this->INTEGRATOR_TOL && Jperr > this->INTEGRATOR_TOL )
    JpInitWarning = true;
  
  // Initialization of eFn
  FTensors pFInv = ttl::inverse( this->pFn );
  this->eFn(i,j) = this->Fn(i,l) * pFInv(l,j);
  this->eFnp1 = this->eFn;
  
  // Assignment of the Second Piola-Kirchoff stress
  this->Sn = S0;
  this->Snp1 = S0;
  
  this->KSn = this->KirchoffStressTensor( this->pFn, this->eFn, this->Sn );
  this->KSnp1 = this->KSn;

  this->sigman = this->CauchyStressTensor( this->eFn, this->Sn );
  this->sigmanp1 = this->sigman;
  
  // Check the compatibility for the Second Piola-Kirchoff stress
  double PKerr = 0, PKmaxerr = 0;
  FTensors SCheck = this->SecondPKTensor( this->eFn, this->pcn );
    
  for ( unsigned q=0; q<dim*dim; q++ )
  {
    PKerr = fabs( SCheck.get(q) - this->Sn.get(q) );
    if (PKerr > PKmaxerr)
      PKmaxerr = PKerr;
  }
  bool PKInitWarning = false;
  if ( PKmaxerr / FrobeniusNorm( this->Sn ) > this->INTEGRATOR_TOL )
    PKInitWarning = true;
  
  // Output
  if ( JpInitWarning )
  {
    std::string Fstr="";
    PrintInVectorForm( this->pFn, Fstr);
    std::cerr << " WARNING: KMS_IJSS2017_Explicit_FE<dim>::Initialize( const double, const FTensors&, const FTensors&, const FTensors&  ) ";
    std::cerr << " Incompatibility between provided pc0 " << this->pcn << " and pF0 " << Fstr << ". \n";
    std::cerr << " Code did not abort but outcomes might be wrong. \n";
  }
  
  if ( PKInitWarning )
  {
    std::string Fstr="", Fpstr="", Sstr="";
    PrintInVectorForm( F0, Fstr);
    PrintInVectorForm( this->pFn, Fpstr);
    PrintInVectorForm( S0, Sstr);
    std::cerr << " WARNING: KMS_IJSS2017_Explicit_FE<dim>::Initialize( const double, const FTensors&, const FTensors&, const FTensors&  ) ";
    std::cerr << " Incompatibility between provided pc0" << pc0 << ", F0 " << Fstr << ", pF0 " << Fpstr << ", and S0 " << Sstr << " according to the eF pF decomposition. \n";
    std::cerr << " Code did not abort but outcomes might be wrong. \n";
  }
  
  
}

template <int dim>
void KMS_IJSS2017_Integration_Algorithms<dim>::set_data_from_PDE(double *Fnp1_in, 
                                                                 double *Fn_in,
                                                                 double *pFnp1_in,
                                                                 double *pFn_in,
                                                                 double pcnp1_in,
                                                                 double pcn_in)
//! Initialize the integrator with pcn, Fn, and Sn at t=0
//! The quantities at n+1 are also initialized pcn, Fn, and Sn at t=0
//! so that the consistent tangent stiffness at the very initial time can be estimated
//! as usual
{

  using namespace ttlindexes;
  memcpy(this->Fnp1.data,  Fnp1_in,  sizeof(double)*dim*dim);
  memcpy(this->Fn.data,    Fn_in,    sizeof(double)*dim*dim);
  memcpy(this->pFnp1.data, pFnp1_in, sizeof(double)*dim*dim);
  memcpy(this->pFn.data,   pFn_in,   sizeof(double)*dim*dim);
      
  // Zeroing
  this->Dpn = ttl::identity(i,j);
  
  // Assignment of pcn
  this->pcn = pcn_in;
  this->pcnp1 = pcnp1_in;
    
  // Initialization of eF
  FTensors pFnI = ttl::inverse(this->pFn);
  this->eFn(i,j) = this->Fn(i,l)*pFnI(l,j);
  
  FTensors pFnp1I = ttl::inverse(this->pFnp1);
  this->eFnp1(i,j) = this->Fnp1(i,l)*pFnp1I(l,j);  
  
  // Assignment of the Second Piola-Kirchoff stress
  this->Sn   = this->SecondPKTensor(this->pFn,   pcn_in);
  this->Snp1 = this->SecondPKTensor(this->pFnp1, pcnp1_in);
  
  this->KSn   = this->KirchoffStressTensor( this->pFn, this->eFn, this->Sn );
  this->KSnp1 = this->KirchoffStressTensor( this->pFnp1, this->eFnp1, this->Snp1);

  this->sigman   = this->CauchyStressTensor( this->eFn, this->Sn);
  this->sigmanp1 = this->CauchyStressTensor( this->eFnp1, this->Snp1);
 
}


template <int dim>
void KMS_IJSS2017_Integration_Algorithms<dim>::set_data_to_PDE(double *pFnp1_in,
                                                               double *pcnp1_in)
//! Initialize the integrator with pcn, Fn, and Sn at t=0
//! The quantities at n+1 are also initialized pcn, Fn, and Sn at t=0
//! so that the consistent tangent stiffness at the very initial time can be estimated
//! as usual
{
  memcpy(pFnp1_in, this->pFnp1.data, sizeof(double)*dim*dim);
  *pcnp1_in = this->pcnp1;
}

template <int dim>
unsigned KMS_IJSS2017_Integration_Algorithms<dim>::FindpcFromJpAtStepnP1( bool Verbose )
//! NR scheme to estimate pc at step (n+1) from Jp at the same step
//! Returns the number of iterations required to convergence, whereas the new outcome is stored in pcnp1
//! The input data is contained in pFnp1 that must be updated before this function is called
{
  
  // Initial checks
  
  double a1=this->Parameters->hr_a1, a2=this->Parameters->hr_a2, l1=this->Parameters->hr_Lambda1, l2=this->Parameters->hr_Lambda2;
  
  double NRTOL = 1e-6;
  double logJp = log( det( this->pFnp1 ) );
  double fpcn = -(  logJp  + a1 * exp( - l1/ this->pcn ) + a2 * exp( - l2/ this->pcn ) );
  
  if ( Verbose && (fpcn < 0) )
  {
    std::cout << " WARNING: KMS_IJSS2017_Explicit<dim>::FindpcFromJp( ) - Unexpected sign for f(pc(n)): it is negative, which means that pc decreased.";
    std::cerr << " Code did not abort but outcomes might be wrong. \n";
  }
  
  double DfDpcn = ( a1 * l1 * exp( - l1/ this->pcn ) + a2 * l2 * exp( - l2/ this->pcn ) ) / ( this->pcn * this->pcn );
  
  // Algorithm selection, either NR or bisections
  
  unsigned it=0;
  
  if ( fabs( DfDpcn ) < NRTOL )
    // the first derivative is almost zero, and the NR scheme becomes very bad conditioned.
    // the zero is found via bisection method
  {
    if ( Verbose )
    {
      std::cout << " WARNING: KMS_IJSS2017_Explicit<dim>::FindpcFromJp( ) - The first derivative is almost zero at pc(n). pc(n+1) will be estimated via bisection method. ";
    }
    
    double a = this->pcn, b = a, fpcnp1 = 1, c=a;
    // seeking for the upper extreme of the bounding interval
    do
    {
      b *= 2;
      fpcnp1 = -(  logJp  + a1 * exp( - l1/ b ) + a2 * exp( - l2/ b ) );
    } while ( fpcnp1 > 0 );
    
    // at this point we have fpcn>0 and fpcnp1 <0
    // the bisection method can run
    unsigned maxIt=500;
    double fpccheck = 1;
    
    while (it < maxIt)
    {
      it++;
      c = 0.5 * (a+b);
      fpccheck = -(  logJp  + a1 * exp( - l1/ c ) + a2 * exp( - l2/ c ) );
      if ( fabs( fpccheck ) < this->INTEGRATOR_TOL )
        break;
      else if ( fpccheck > 0 )
        a = c;
      else
        b =c;
    }
    
    this->pcnp1 = c;
    
    if ( Verbose && it == maxIt )
    {
      std::cerr << " WARNING: KMS_IJSS2017_Explicit<dim>::FindpcFromJp( ) - Convergence was not achieved in the bisection scheme after " << maxIt << " iterations.";
      std::cerr << " Code did not abort but outcomes might be wrong. \n";
    }
  }
  
  else
    // the NR method is used to find pc(n+1)
  {
    unsigned maxIt=100;
    this->pcnp1 = this->pcn;
    
    while (it < maxIt) {
      
      it++;
      double dpc = -(  logJp  + a1 * exp( - l1/ this->pcnp1 ) + a2 * exp( - l2/ this->pcnp1 ) ) * ( this->pcnp1 * this->pcnp1 );
      dpc /= a1 * l1 * exp( - l1/ this->pcnp1 ) + a2 * l2 * exp( - l2/ this->pcnp1 ) ;
      this->pcnp1 += dpc;
      
      if ( fabs( logJp  + a1 * exp( - l1/ this->pcnp1 ) + a2 * exp( - l2/ this->pcnp1 ) ) < this->INTEGRATOR_TOL )
        break;
    }
    
    if ( Verbose && it == maxIt )
    {
      std::cerr << " WARNING: KMS_IJSS2017_Explicit<dim>::FindpcFromJp( ) - Convergence was not achieved in the Newton-Raphson scheme after " << maxIt << " iterations.";
      std::cerr << " Code did not abort but outcomes might be wrong. \n";
    }
  }
  
  return it;
}



template <int dim>
unsigned KMS_IJSS2017_Integration_Algorithms<dim>::SecondPKTensorAtStepnP1( bool Verbose )
//! Estimate of the second Piola-Kirchoff stress from data at step (n+1)
//! The volumetric part is peculiar of this model, whereas the isochoric contribution is the usual Neo-Hookean
//! with material parameters that are functions of the internal variable pc
{
  
  
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'l'> l;
  
  FTensors Id2 = ttl::identity(i,j); // ttl::Delta<2,dim,double>(1);

  // Je at step n+1
  double eJnp1 = det( this->eFnp1 ) ;
  
  // Left Cauchy-Green tensor and its inverse at step n+1
  FTensors Ce;
  Ce(i,j) = this->eFnp1(l,i) * this->eFnp1(l,j);
  
  FTensors InvCe = ttl::inverse( Ce );
  
  if ( Verbose )
  {
    std::string Fstr="";
    PrintInVectorForm( InvCe , Fstr );
    std::cout << " InvCe(n+1)=" << Fstr;
  }
  
  // Volumetric contribution
  double exponent = - pow( ( 1.0 - this->c(this->pcnp1) / this->Parameters->c_inf  ), this->Parameters->pl_n  ) / this->Parameters->K_kappa ;
  this->Snp1(i,j) = (
                     0.5 * this->bulkmodulus( this->pcnp1 ) * ( eJnp1 * log(eJnp1) + eJnp1 - 1.0 )
                     +  this->c( this->pcnp1 )
                     - ( this->Parameters->K_p0 + this->c( this->pcnp1 ) ) * pow( eJnp1, exponent )
                     )
  * InvCe(i,j) ;
  
  if ( Verbose )
  {
    std::cout << " this->bulkmodulus( this->pcnp1 )=" << std::setw(20) << std::setprecision(15) << this->bulkmodulus( this->pcnp1 ) << "\n";
    std::cout << " ( eJnp1 * log(eJnp1) + eJnp1 - 1.0 )=" << std::setw(20) << std::setprecision(15) << ( eJnp1 * log(eJnp1) + eJnp1 - 1.0 ) << "\n";
    std::cout << " this->c( this->pcnp1 )=" << std::setw(20) << std::setprecision(15) << this->c( this->pcnp1 ) << "\n";
    std::cout << " ( this->Parameters->K_p0 + this->c( this->pcnp1 ) )=" << std::setw(20) << std::setprecision(15) << ( this->Parameters->K_p0 + this->c( this->pcnp1 ) ) << "\n";
    std::cout << " exponent=" << std::setw(20) << std::setprecision(15) << exponent << "\n";
    std::cout << " pow( eJnp1, exponent )=" << std::setw(20) << std::setprecision(15) << pow( eJnp1, exponent ) << "\n";
    std::string Fstr="";
    PrintInVectorForm( this->Snp1 , Fstr );
    std::cout << " Volumetric Second Piola-Kirchoff S(n+1)=" << Fstr;
  }
  
  // Isochoric contribution
  this->Snp1(i,j) += this->shearmodulus( this->pcnp1 ) * pow( eJnp1, -2.0/3.0 ) * ( Id2(i,j) - Trace(Ce) / 3.0 * InvCe(i,j) );
  
  return 0;
}



template <int dim>
void KMS_IJSS2017_Integration_Algorithms<dim>::IntegratorTest ( std::function<void(const FTensors&, double, FTensors&)> HistoryOfF, bool Verbose )
//! This method tests the KMS_IJSS2017_Explicit_FE integrator
//! 0 - Set the time frame
//! 1 - Initial conditions
//! 2 - Update until the final time is overcome
{
  
  std::string str="";
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  
  this->IO->log().precision(3);
  
  // 0
  // Set the time-frame
  // this steps has been replaced by the constructor and the definition of the Time Integration Data Manager
  // this->SetTimeFrame( 0.0, FinalTime, Dt, 0.0, 0 );
  
  // 1
  // Initialization parameters
  double p0 = this->Parameters->K_p0 ;
  double HardLawJp0Coeff = pow(exp( this->HardeningLaw(p0) ), 1.0/(double) dim);
  
  FTensors Id2 = ttl::identity(i,j);
  FTensors F0 = Id2(i,j) * HardLawJp0Coeff ;
  FTensors pF0 = F0;
  FTensors eF0 = Id2(i,j);
  
  FTensors S0 = this->SecondPKTensor( eF0, p0 );
  
  this->Initialize( p0, F0, pF0, S0 );
  
  if ( Verbose )
  {
    this->AsAString( str, Verbose );
    std::cout << "\n  KMS_IJSS2017 Integrator Initialization\n";
  }
  
  this->IO->log() << "\n  Time = " << std::fixed << std::setprecision(3) << this->TimeIntegrationData->CurrentTime() << ", step = " << this->TimeIntegrationData->TimeStep() << "\n" << std::flush;
  
  
  // 2
  // this->CurrentTime = InitialTime;
  //unsigned step = 0;  - replaced by the definition of timestep in the base class

  FTensors F;
  
  this->PrintStep( this->TimeIntegrationData->TimeStep(), this->TimeIntegrationData->CurrentTime() );
  do {
    
    // update the step and check if the time overcomes the
    // final time
    
    if ( this->TimeIntegrationData->NewStep() > 0 )
    {
      
      if ( Verbose )
        std::cout << "\n  Time = " << std::fixed << std::setprecision(2) << this->TimeIntegrationData->CurrentTime()  << ", step = " << std::setw(8) << std::setprecision(2) << this->TimeIntegrationData->TimeStep()  << "\n";
      if ( this->IO->StepToPrint( this->TimeIntegrationData->TimeStep()  ) )
        this->IO->log() << "  Time = " << std::fixed << std::setprecision(3) << this->TimeIntegrationData->CurrentTime()  << ", step = "  << this->TimeIntegrationData->TimeStep() << " running" << std::flush;

      HistoryOfF( F0, this->TimeIntegrationData->CurrentTime(), F );

      if ( Verbose )
        this->VerboseStepUpdate( F, this->TimeIntegrationData->Deltat() , Verbose );
      else
        this->StepUpdate( F, this->TimeIntegrationData->Deltat() , Verbose );
      
      this->StepClose();
      this->PrintStep( this->TimeIntegrationData->TimeStep(), this->TimeIntegrationData->CurrentTime()  );

      if ( this->IO->StepToPrint( this->TimeIntegrationData->TimeStep() ) )
        this->IO->log() << " / completed. \n" << std::flush;

    }
    else
      break;
    
  } while ( true );

  
}





template <int dim>
void KMS_IJSS2017_Integration_Algorithms<dim>::ConsistentTangentTest ( std::function<void(const FTensors&, double, FTensors&)> HistoryOfF, bool Verbose )
//! This method tests the KMS_IJSS2017_Explicit_FE integrator
//! 0 - Set the time frame
//! 1 - Initial conditions
//! 2 - Update until the final time is overcome
{
  
  using namespace ttlindexes;
  
  std::string str="";
  
  this->IO->log().precision(3);
  
  // 0
  // Set the time-frame
  // this steps has been replaced by the constructor and the definition of the Time Integration Data Manager
  // this->SetTimeFrame( 0.0, FinalTime, Dt, 0.0, 0 );
  
  // 1
  // Initialization parameters
  double p0 = this->Parameters->K_p0 ;
  double HardLawJp0Coeff = pow(exp( this->HardeningLaw(p0) ), 1.0/(double) dim);
  
  FTensors Id2 = ttl::identity(i,j);
  FTensors F0 = Id2(i,j) * HardLawJp0Coeff ;
  FTensors pF0 = F0;
  FTensors eF0 = Id2(i,j);
  
  FTensors S0 = this->SecondPKTensor( eF0, p0 );
  
  this->Initialize( p0, F0, pF0, S0 );
  
  if ( Verbose )
  {
    this->AsAString( str, Verbose );
    std::cout << "\n  KMS_IJSS2017 Integrator Initialization\n";
    
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

    PrintInVectorForm( this->KSnp1 , Fstr );
    std::cout << "  Kirchoff KSnp1 = " << Fstr;
    Fstr="";

    PrintInVectorForm( this->sigmanp1 , Fstr );
    std::cout << "  Cauchy sigmanp1 = " << Fstr;
    Fstr="";

  }
  
  this->IO->log() << "\n  Time = " << std::fixed << std::setprecision(3) << this->TimeIntegrationData->CurrentTime() << ", step = " << this->TimeIntegrationData->TimeStep() << "\n" << std::flush;
  
  
  // 2
  FTensors F;
  
  this->PrintStep( this->TimeIntegrationData->TimeStep(), this->TimeIntegrationData->CurrentTime() );
  do {
    
    // update the step and check if the time overcomes the
    // final time
    
    if ( this->TimeIntegrationData->NewStep() > 0 )
    {
      
      if ( Verbose )
        std::cout << "\n  Time = " << std::fixed << std::setprecision(2) << this->TimeIntegrationData->CurrentTime()  << ", step = " << std::setw(8) << std::setprecision(2) << this->TimeIntegrationData->TimeStep()  << "\n";
      if ( this->IO->StepToPrint( this->TimeIntegrationData->TimeStep() ) )
        this->IO->log() << "  Time = " << std::fixed << std::setprecision(3) << this->TimeIntegrationData->CurrentTime()  << ", step = "  << this->TimeIntegrationData->TimeStep() << " running" << std::flush;
      
      HistoryOfF( F0, this->TimeIntegrationData->CurrentTime(), F );
      
      // consistent tangent test
      ttl::Tensor<4, dim, double> dmdf = {};
      ttl::Tensor<4, dim, double> dsdf = {};
      this->DMDFandDSDF( dmdf, dsdf, this->TimeIntegrationData->Deltat() , Verbose );
      
      if ( Verbose )
      {
        std::string Fstr="";
        
        PrintInVectorForm( dmdf , Fstr );
        std::cout << "  dmdf = " << Fstr << "\n";
        Fstr="";
        
      }
      
      // step update
      if ( Verbose )
        this->VerboseStepUpdate( F, this->TimeIntegrationData->Deltat() , Verbose );
      else
        this->StepUpdate( F, this->TimeIntegrationData->Deltat() , Verbose );
      
      this->StepClose();
      this->PrintStep( this->TimeIntegrationData->TimeStep(), this->TimeIntegrationData->CurrentTime()  );
      
      if ( this->IO->StepToPrint( this->TimeIntegrationData->TimeStep() ) )
        this->IO->log() << " / completed. \n" << std::flush;
      
    }
    else
      break;
  } while ( true );

  
  
}





#endif /* KMS_IJSS2017_integrator_templates_h */
