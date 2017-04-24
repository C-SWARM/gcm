//
//  KMS-IJSS2017-templates.hpp
//  ttl-learning
//
//  Created by Alberto Salvadori on 1/5/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

/*!
 
 \file KMS-IJSS2017-explicit-integrator-templates.h
 \brief Template code file for the class base of explicit integration algorithms for the KMS-IJSS2017 model
 
 This file contains template classes and methods that pertain to every explicit integrator of the poroviscoplastic model described in the paper
 "A. Krairi, K. Matouš, A. Salvadori, A poro-viscoplastic constitutive model for granular materials at finite strain, submitted to IJSS (2017)".
 
 */


#ifndef KMS_IJSS2017_Explicit_FE_templates_hpp
#define KMS_IJSS2017_Explicit_FE_templates_hpp


// include

#include <stdio.h>
#include <sstream>
#include <random>
#include <iostream>
#include <iomanip>
#include <math.h>



// class KMS_IJSS2017_Explicit_FE
// forward euler integrator
// ******************************


template <int dim>
KMS_IJSS2017_Explicit_FE<dim>::KMS_IJSS2017_Explicit_FE(KMS_IJSS2017_Parameters* P, KMS_IJSS2017_IO* InpOut, TimeIntegrationDataManager* TIDM )
: KMS_IJSS2017_Integration_Algorithms<dim>( P, InpOut, TIDM )
//! Constructor with material parameters and IO
{}


// Public methods
// --------------


template <int dim>
void KMS_IJSS2017_Explicit_FE<dim>::AsAString( std::string& str, bool Verbose )
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
    str += "\n\n";
    str += "KMS_IJSS2017 Explicit Integrator for dimension ";
    str += std::to_string(dim);
    str += ". \n";
  }
  
  // data output in vector form
  str += " Deltat = "; str += std::to_string(this->Deltat); str += "\n";
  str += " pc(n) = "; str += std::to_string(this->pcn);
  str += ", pc(n+1) = "; str += std::to_string(this->pcnp1); str += "\n";
  str += " F(n) = "; PrintInVectorForm( this->Fn, str );
  str += " F(n+1) = "; PrintInVectorForm( this->Fnp1, str );
  str += " Fe(n) = "; PrintInVectorForm( this->eFn, str );
  str += " Fe(n+1) = "; PrintInVectorForm( this->eFnp1, str );
  str += " Fp(n) = "; PrintInVectorForm( this->pFn, str );
  str += " Fp(n+1) = "; PrintInVectorForm( this->pFnp1, str );
  str += " S(n) = "; PrintInVectorForm( this->Sn, str );
  str += " S(n+1) = "; PrintInVectorForm( this->Snp1, str );
  str += " Dp(n) = "; PrintInVectorForm( this->Dpn, str );
  
}





template <int dim>
void KMS_IJSS2017_Explicit_FE<dim>::TestSuite( std::string& str )
//! This method tests some of the features of the KMS_IJSS2017_Explicit_FE integrator
{
  
  // Inverse method test
  
  double lower_bound = 0;
  double upper_bound = 10000;
  std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
  
  static constexpr ttl::Index<'i'> i;
  static constexpr ttl::Index<'j'> j;
  static constexpr ttl::Index<'l'> l;

  double maxerr = 0;
  int maxlp = 100;
  
  for (unsigned lp=0; lp<maxlp; lp++)
  {
    std::default_random_engine re( (unsigned) std::time(0) );
    
    FTensors A;
    for (unsigned ii=0; ii<dim; ii++)
      for (unsigned jj=0; jj<dim; jj++)
        A[i][j] = unif(re);
    
    FTensors InvA = ttl::inverse(A);
    //Inv(A, InvA);
    
    FTensors Id;
    Id(i,j) = A(i,l) * InvA(l,j);
    double err = fabs( FrobeniusNorm( Id )/sqrt(dim) - 1 );
    
    if (err > maxerr) maxerr = err;
  }
  
  std::ostringstream outs;
  outs << " Inverse method test for dimension " << dim << ". Max error " << std::setw(20) << std::setprecision(15)  << maxerr  << " on the Frobenius norm on " << maxlp << " attempts \n";
  
  
  // FindpcFromJp method test
  
  unsigned it=0;
  this->pcn = 30;
  double pc0 =  this->pcn;
  double x = pow( exp(-0.4), 1.0/dim );
  
  ttl::Tensor<2, dim, double> Id2 = ttl::identity( i,j ); // ttl::Delta<2,dim,double>(1);
  this->pFnp1 = Id2(i,j) * x; // ttl::Delta<2,dim,double>(x); //

  it = this->FindpcFromJpAtStepnP1();
  
  outs << " FindpcFromJp method test for dimension " << dim << ". log(jp) = " << std::setw(20) << std::setprecision(15) << log(det(this->pFnp1 )) << ", pcn= " << std::setw(20) << std::setprecision(15)  << pc0
       << ", pcnp1= " << std::setw(20) << std::setprecision(15) << this->pcnp1
       << ", TOL= " << this->TOL << ", it= " << it << "\n";

  str +=  outs.str();
  
  
}








template <int dim>
unsigned KMS_IJSS2017_Explicit_FE<dim>::StepUpdate( const FTensors& updF, const double dt, bool Verbose )
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
  
  if ( Verbose )
  {
    std::string Fstr="";
    PrintInVectorForm( Kirchoffn , Fstr );
    std::cout << " Kirchoff (n)=" << Fstr;
    Fstr="";
    PrintInVectorForm( DevKirchoffn , Fstr );
    std::cout << " Dev( Kirchoff (n) )=" << Fstr;
    std::cout << " Tr( Kirchoff (n) )=" << std::setw(20) << std::setprecision(15)  << pi << "\n";
  }

  
  // 1
  this->Fnp1 = updF;
  
  if ( Verbose )
  {
    std::string Fstr="";
    PrintInVectorForm( this->Fnp1 , Fstr );
    std::cout << " F(n+1)=" << Fstr;
  }
  
  // 2
  double gmdv = this->gammadot_v(pi, this->pcn);
  double gmdd = this->gammadot_d(FrobNormDevKirchoffn, this->pcn);
  
  if ( std::fpclassify(gmdv) == FP_ZERO && std::fpclassify(gmdd) == FP_ZERO)
  // if ( fabs(gmdv) < this->ZERO && fabs(gmdd) < this->ZERO )
  {
    if ( Verbose )
    {
      std::cout << " pc(n) = " << std::setw(20) << std::setprecision(15)  << this->pcn  << ", pi_m(pc) = " << std::setw(20) << std::setprecision(15)  << this->pi_m(this->pcn)
                << ", gammadot_v = " << std::setw(20) << std::setprecision(15) << gmdv << ", gammadot_d= " << std::setw(20) << std::setprecision(15) << gmdd << " Elastic step. \n";
    }
    
    this->Dpn = Id2(i,j) * 0 ; // ttl::Delta<2,dim,double>(0);
    this->pFnp1(i,j) =  this->pFn(i,j);
    this->pcnp1 = this->pcn;
    
  }
  else
  {
  
    if ( Verbose )
    {
      std::cout << " pc(n) = " << std::setw(20) << std::setprecision(15) << this->pcn
                << ", pi_m(pc) = "  << std::setw(20) << std::setprecision(15)  << this->pi_m(this->pcn)
                << ", gammadot_v = " << std::setw(20) << std::setprecision(15)  << gmdv
                << ", gammadot_d= " << std::setw(20) << std::setprecision(15) << gmdd << "\n";
    }
    
    // 3
    FTensors Nn = dev( this->Sn );
    double NnNorm= FrobeniusNorm( Nn );
    
    if ( NnNorm > this->INTEGRATOR_TOL )
      this->Dpn(i,j) = gmdd / NnNorm * Nn(i,j) + 1.0/3.0 * ( this->betaD( this->pcn ) * gmdd - this->betaC( this->pcn ) * gmdv ) * Id2(i,j);
    else
      this->Dpn(i,j) = 1.0/3.0 * ( this->betaD( this->pcn ) * gmdd - this->betaC( this->pcn ) * gmdv ) * Id2(i,j);
    
    if ( Verbose )
    {
      std::string Fstr="";
      PrintInVectorForm( this->Dpn , Fstr );
      std::cout << " Plastic strectch Dp(n+1)=" << Fstr;
      std::cout << " Tr(Dp(n+1)) = " << std::setw(20) << std::setprecision(15)  << Trace(this->Dpn) << "\n";
    }
    
    // 4
    this->pFnp1(i,j) =  this->pFn(i,j) + this->Dpn(i,l) * this->pFn(l,j) * this->TimeIntegrationData->Deltat();
    
    if ( Verbose )
    {
      std::string Fstr="";
      PrintInVectorForm( this->pFnp1 , Fstr );
      std::cout << " Plastic strain Fp(n+1)=" << Fstr;
      Fstr="";
      std::cout << " Jp = " << std::setw(20) << std::setprecision(15)  << det( this->pFnp1 ) << ", log(Jp) = " << std::setw(20) << std::setprecision(15) << log( det( this->pFnp1 ) ) << "\n";
    }
    
    // 5
    this->FindpcFromJpAtStepnP1( Verbose );
    
    if ( Verbose )
    {
      std::cout << " pc(n+1)=" << std::setw(20) << std::setprecision(15) << this->pcnp1 << "\n";
    }
    
  }

  // 6
  FTensors buffer = ttl::inverse( this->pFnp1 );
  // Inv( this->pFnp1, buffer );
  this->eFnp1(i,j) =  this->Fnp1(i,l) * buffer(l,j);
  
  if ( Verbose )
  {
    std::string Fstr="";
    PrintInVectorForm( this->eFnp1 , Fstr );
    std::cout << " Elastic strain Fe(n+1)=" << Fstr;
    std::cout << " Je(n+1)=" << det( this->eFnp1 ) << "\n";
    
    // Left Cauchy-Green tensor and its inverse at step n+1
    FTensors Ce;
    Ce(i,j) = this->eFnp1(l,i) * this->eFnp1(l,j);
    Fstr="";
    PrintInVectorForm( Ce , Fstr );
    std::cout << " Left Cauchy-Green tensor Ce(n+1)=" << Fstr;
  }

  // 7
  this->SecondPKTensorAtStepnP1( Verbose );
  
  if ( Verbose )
  {
    std::string Fstr="";
    PrintInVectorForm( this->Snp1 , Fstr );
    std::cout << " Second Piola-Kirchoff S(n+1)=" << Fstr << "\n";
  }

  // 8
  this->KSnp1 = this->KirchoffStressTensor( this->pFnp1, this->eFnp1, this->Snp1 );
  
  if ( Verbose )
  {
    std::string Fstr="";
    PrintInVectorForm( this->KSnp1 , Fstr );
    std::cout << " Kirchoff stress S(n+1)=" << Fstr << "\n";
  }


  // 9
  this->sigmanp1 = this->CauchyStressTensor( this->eFnp1, this->Snp1 );
  
  if ( Verbose )
  {
    std::string Fstr="";
    PrintInVectorForm( this->sigmanp1 , Fstr );
    std::cout << " Cacuhy stress S(n+1)=" << Fstr << "\n";
  }

  return 1;
  
}











#endif /* KMS_IJSS2017_Explicit_FE_templates_hpp */
