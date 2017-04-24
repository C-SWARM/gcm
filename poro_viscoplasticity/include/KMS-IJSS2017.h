//
//  KMS-IJSS2017.h
//  ttl-learning
//
//  Created by alberto salvadori on 12/13/16.
//  Copyright © 2016 alberto salvadori. All rights reserved.
//

/*!
 
 \file KMS-IJSS2017.h
 \brief Header file for the KMS-IJSS2017 model
 
 This file contains classes and methods for the poroviscoplastic model desccribed in the paper
 "A. Krairi, K. Matouš, A. Salvadori, A poro-viscoplastic constitutive model for granular materials at finite strain, submitted to IJSS (2017)".
 
 */


#ifndef KMS_IJSS2017_h
#define KMS_IJSS2017_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ttl/ttl.h>

#include "MechanicalModels.h"
#include "FeFp_integrator.h"
#include "FeFpModels.h"
#include "TimeIntegration_Manager.h"


class KMS_IJSS2017_IO : public FeFpModels_IO
/*!
 \brief This class contains and handles the input and output for the KMS-IJSS2017 model
 */
{
protected:
  
  std::ofstream pcfile;     //!< output of the inner variable pc for KMS_IJSS2017
  std::ofstream matfile;    //!< output of the material properties for KMS_IJSS2017
  
  //unsigned PrintEveryNSteps;
  
  
public:
  // Constructors & destructors
  KMS_IJSS2017_IO(){}
  KMS_IJSS2017_IO( const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, unsigned = 1 );
  virtual ~KMS_IJSS2017_IO();
  
  // Accessing IO
  std::ofstream& pc(){ return pcfile; }
  std::ofstream& mat(){ return matfile; }
  
};



class KMS_IJSS2017_Parameters
/*!
 \brief This class contains and handles the parameters for the KMS-IJSS2017 model
 */
{
  template <int dim> friend class KMS_IJSS2017;
  template <int dim> friend class KMS_IJSS2017_Integration_Algorithms;
  template <int dim> friend class KMS_IJSS2017_Explicit_FE;
  template <int dim> friend class KMS_IJSS2017_Implicit;
  template <int dim> friend class KMS_IJSS2017_Implicit_BE;
  template <int dim> friend class KMS_IJSS2017_Implicit_BE_Staggered;

private:
  
  // Material parameters
  double yf_M, yf_alpha;                          //!< 1. Yield function parameters
  double flr_m, flr_gamma0;                       //!< 2. Flow rule parameters
  double hr_a1, hr_a2, hr_Lambda1, hr_Lambda2;    //!< 3. Hardening rule parameters
  double c_inf, c_Gamma;                          //!< 4. Cohesion rule parameters
  double d_B, d_pcb;                              //!< 5. Transition rule parameters
  double mu_0, mu_1;                              //!< 6. Shear modulus parameters
  double K_p0, K_kappa;                           //!< 7. Bulk modulus parameters
  double pl_n;                                    //!< 8. Power law exponent
  double cf_g0, cf_pcinf;                         //!< 9. Compaction function parameters
  
  bool smMbrackets;                               //!< smooth Mackauley brackets flag
  
public:
  
  // Constructors
  KMS_IJSS2017_Parameters                         //!< Constructor with material parameters
          (double, double, double, double, double, double, double, double,
           double, double, double, double, double, double, double, double,
           double, double, double, bool);
  
  // Destructors
  virtual ~KMS_IJSS2017_Parameters(){}
  
  // Methods
  void Checks( bool );                            //!< This method performs some checks and prints warnings,
  
  // IO
  void AsAString( std::string& );                 //!< This method prints the model features as a string
  
};





template <int dim>
class KMS_IJSS2017_Integrator : public FeFp_Integrator<dim>
/*!
 \brief This abstract class contains data and methods for the integration of the KMS-IJSS2017 model
 */
{
  
  typedef ttl::Tensor<2, dim, double> FTensors;
  
protected:

  // IO
  KMS_IJSS2017_IO* IO;                              //!< Input-output manager. It is created elsewhere and passed as a pointer
  bool Print;                                       //!< Flag for printing the solution at each step on files.
  void PrintData( size_t, double, bool = false );   //!< Prints data of each step on the proper files

  // Internal variables
  double pcn, pcnp1;                                //!< Internal variable at time steps n and n+1
  
  // pc == pc_inf flag
  bool pcEQpc_inf;

  // Methods
  void StepClose( );                                //!< Assigns data at step n from data at step n+1
  

public:
  
  // Constructors
  KMS_IJSS2017_Integrator();                         //!< Default Constructor
  KMS_IJSS2017_Integrator( KMS_IJSS2017_IO* , TimeIntegrationDataManager*);       //!< Constructor with IO
  

};


template <int dim>
class KMS_IJSS2017 : public MechanicalModels<dim>
/*!
 \brief This class contains data and methods for the KMS-IJSS2017 model
 */
{
  
  typedef ttl::Tensor<2, dim, double> FTensors;

protected:
  KMS_IJSS2017_Parameters* Parameters;              //!< Material parameters. They are created elsewhere and passed as a pointer
  
  // Methods
  double a(double);                                 //!< Yield function semi-axis
  double b(double, double);                         //!< Yield function ellipses semi-axis as a function of pc and of pressure, pi
  double betaC(double);                             //!< Compaction function
  double betaD(double);                             //!< Shear-induced dilatancy function
  double c(double);                                 //!< Cohesion
  double d(double);                                 //!< Transition function
  double gammadot_d(double,double);                 //!< Equivalent shear plastic strain rate as a function of tau and pc
  double gammadot_v(double,double);                 //!< Equivalent volumetric plastic strain rate as a function of pi and pc
  double g_pi(double, double);                      //!< "Yield-like" pressure point
  double g_tau(double);                             //!< "Yield-like" shear point
  double pi_m(double);                              //!< Yield function ellipses centroid as a function of pc
  
  double shearmodulus(double);                      //!< Shear modulus mu as a function of pc
  double bulkmodulus(double);                       //!< Bulk modulus K as a function of pc

  double DcDpc(double);                             //!< Derivative of the cohesion wrt pc
  double DdDpc(double);                             //!< Derivative of the transition function wrt pc
  double DmuDpc(double);                            //!< Derivative of the Shear modulus wrt pc
  double DbulkDpc(double);                          //!< Derivative of the Bulk modulus wrt pc
  double DalphaDpc(double);                         //!< Derivative of the coefficient alpha wrt pc
  double Dg_tauDpc(double);                         //!< Derivative of the coefficient g_tau wrt pc
  double Dpi_mDpc(double);                          //!< Derivative of the coefficient pi_m wrt pc
  double Dg_piDpc(double, double);                  //!< Derivative of the coefficient g_pi wrt pc

public:
  
  // Constructors
  KMS_IJSS2017()  { Parameters=NULL;}               //!< Default Constructor
  KMS_IJSS2017( KMS_IJSS2017_Parameters* );          //!< Constructor with material parameters
  
  // Destructors
  virtual ~KMS_IJSS2017(){}
  
  // Generic Methods
  double HardeningLaw( double );                                                              //!< Hardening law, returns Log(jp) as a function of pc
  double DHardeningLawDpc( double );                                                          //!< Derivative of the Hardening law with respect to pc
  FTensors SecondPKTensor( const FTensors&, const double );                                   //!< Estimate of the second Piola-Kirchoff stress from a given eF and a given pc
  FTensors KirchoffStressTensor( const FTensors&, const FTensors&, const FTensors& );         //!< Estimate of the Kirchoff stress from a given pF, eF, S
  FTensors CauchyStressTensor( const FTensors&, const FTensors& );                            //!< Estimate of the Cauchy stress from a given eF, S
  
  
};



template <int dim>
class KMS_IJSS2017_Integration_Algorithms  : public KMS_IJSS2017<dim>, public KMS_IJSS2017_Integrator<dim>
/*!
 \brief This class containes data and methods for the KMS-IJSS2017 model and the algorithms shared by all integrators
 */
{
  typedef ttl::Tensor<2, dim, double> FTensors;

  
protected:
  
  // IO
  void PrintStep( size_t, double, bool = false );               //!< Prints data and material properties of each step on the proper files

public:
  
  // Constructors
  KMS_IJSS2017_Integration_Algorithms() : KMS_IJSS2017<dim> () {}                       //!< Default Constructor
  KMS_IJSS2017_Integration_Algorithms(KMS_IJSS2017_Parameters* , KMS_IJSS2017_IO*, TimeIntegrationDataManager*);      //!< Constructor with material parameters and IO
  
  // Destructors
  virtual ~KMS_IJSS2017_Integration_Algorithms(){}
  
  // Methods shared with all integrators
  void Initialize( const double, const FTensors&, const FTensors&,  const FTensors& );  //!< Initialize the integrator with pcn, Fn, pFn, and Sn at t=0
  unsigned FindpcFromJpAtStepnP1( bool = false);                 //!< Estimate pc at step (n+1) from Jp at the same step, (n+1)
  unsigned SecondPKTensorAtStepnP1( bool = false );              //!< Estimate the second Piola-Kirchoff stress at step (n+1) from data at the same time step, (n+1)

  // Integrator methods
  virtual unsigned StepUpdate( const FTensors&, const double, bool ){ return 0; }             //!< Updates stresses and internal variables from a given FnP1 and dt
  virtual unsigned VerboseStepUpdate( const FTensors&, const double, bool ){ return 0; }      //!< Updates stresses and internal variables from a given FnP1 and dt

  // Integrator methods: iterative DMDF used for debug
  virtual ttl::Tensor<4, dim, double> IterativeDMDF(ttl::Tensor<4, dim, double>, const double, const bool) { using namespace ttlindexes; return ttl::zero(i,j,k,l); }                              //!< This function evaluates iteratively the fourth order tensor DM/DF
  virtual ttl::Tensor<4, dim, double> DMDF(const double, const bool){ using namespace ttlindexes; return ttl::zero(i,j,k,l); }                                                                     //!< This function evaluates the fourth order tensor DM/DF
  virtual void DMDFandDSDF(ttl::Tensor<4, dim, double> &, ttl::Tensor<4, dim, double> &, const double, const bool){}                      //!< This function evaluates the fourth order tensors DM/DF and DS/DF (total derivatives)
  

  // Test methods
  void IntegratorTest( std::function<void(const FTensors&, double, FTensors&)>, bool Verbose = false );  //!< This method tests the KMS_IJSS2017_Integration_Algorithms integrator
  void ConsistentTangentTest( std::function<void(const FTensors&, double, FTensors&)>, bool Verbose = false );  //!< This method tests the KMS_IJSS2017_Integration_Algorithms integrator
};



template <int dim>
class KMS_IJSS2017_Explicit_FE  : public KMS_IJSS2017_Integration_Algorithms<dim>
/*!
 \brief This class containes data and methods for the KMS-IJSS2017 model and the algorithms for the explicit forward euler integrator
 */
{
  
  typedef ttl::Tensor<2, dim, double> FTensors;

public:

  // Constructors
  KMS_IJSS2017_Explicit_FE() : KMS_IJSS2017_Integration_Algorithms<dim> () {}                       //!< Default Constructor
  KMS_IJSS2017_Explicit_FE(KMS_IJSS2017_Parameters* , KMS_IJSS2017_IO*, TimeIntegrationDataManager*);      //!< Constructor with material parameters and IO
  
  // Destructors
  virtual ~KMS_IJSS2017_Explicit_FE(){}
  
  // IO
  void AsAString( std::string&, bool );

  // Integrator methods
  unsigned StepUpdate( const FTensors&, const double, bool );      //!< Updates stresses and internal variables from a given FnP1 and dt
  
  // Test methods
  // void IntegratorTest( std::function<void(const FTensors&, double, FTensors&)>,
  //                     double, double, double, bool Verbose = false );  //!< This method tests the KMS_IJSS2017_Explicit_FE integrator
  void TestSuite( std::string& );                                      //!< This method tests some of the features of the KMS_IJSS2017_Explicit_FE integrator

};




template <int dim>
class KMS_IJSS2017_Implicit_BE : public KMS_IJSS2017_Integration_Algorithms<dim>
/*!
 \brief This class contains data and methods for the KMS-IJSS2017 model and the algorithms for the Backward Euler implicit integrator
 */
{
  
  typedef ttl::Tensor<2, dim, double> FTensors;
  
protected:
  
  FTensors UpdateZwx( FTensors );             //!< Updates operator Zwx - see page 17b of the Consistent tangent operator notes
  FTensors UpdateHattedZwx( FTensors );       //!< Updates operator Zwx - see page 17b of the Consistent tangent operator notes, See page 1 of the "Debugging and code optimization" notes
  FTensors DSDpc();                           //!< Derivative of S wrt pc - see page 23 of the Consistent tangent operator notes
  ttl::Tensor<4, dim, double> dSdF( const FTensors&, const FTensors&, const FTensors&, const double );         //!< Derivative of S wrt F (partial derivative) - see page 25 of the Consistent tangent operator notes
  ttl::Tensor<4, dim, double> Minor_cdkl( const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, ttl::Tensor<4, dim, double> );                      //!< Fourth order tensor Minor_cdkl() - see page 17 of the Consistent tangent operator notes
  ttl::Tensor<4, dim, double> Major_cdwx( const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, ttl::Tensor<4, dim, double> );                      //!< Fourth order tensor Major_cdwx() - see page 17b of the Consistent tangent operator notes
  ttl::Tensor<4, dim, double> HattedMajor_cdwx( const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, ttl::Tensor<4, dim, double>, const double );  //!< Fourth order tensor HattedMajor_cdwx() - see page 17b of the Consistent tangent operator notes and page 2 of the "Debugging aand code optimization" notes
  ttl::Tensor<4, dim, double> Phi_abzy( const FTensors&, bool );                      //!< Fourth order tensor Phi_abzy() - see page 18 of the Consistent tangent operator notes
  
  FTensors DgammadotdDM(const ttl::Tensor<4, dim, double>&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&,  const FTensors&, const double );    //!< Derivative of gammadotd wrt to M - see notes
  FTensors DgammadotdDM(const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const double );    //!< Derivative of gammadotd wrt to M - see notes
  FTensors DnormdevtauDtau( const FTensors&, const bool = false );                                                              //!< Derivative of the Frobenus norm of dev(A) wrt to A
  ttl::Tensor<4, dim, double> DKSDM( const ttl::Tensor<4, dim, double>&, const FTensors&, const FTensors& );                    //!< Derivative of the Kirchoff Stress Tensor wrt to M - see notes
  ttl::Tensor<4, dim, double> DGammaDM( const ttl::Tensor<4, dim, double>&, const FTensors&, const FTensors& );                 //!< Derivative of the Tensor Gamma wrt to M - see notes
  ttl::Tensor<4, dim, double> HattedDGammaDM( const ttl::Tensor<4, dim, double>&, const double, const FTensors&, const FTensors& );           //!< Derivative of the Tensor Gamma wrt to M - see notes and and page 2 of the "Debugging aand code optimization" notes
  ttl::Tensor<4, dim, double> DSDM( const FTensors&, const FTensors&, const FTensors&, const double );                          //!< Derivative of the Tensor S (Second Piola-Kirchoff stress) wrt to M - see notes
  ttl::Tensor<4, dim, double> DSisoDM( const FTensors&, const FTensors&, const FTensors&, const double );                       //!< Derivative of the isochoric part of S (Second Piola-Kirchoff stress) wrt to M - see notes
  ttl::Tensor<4, dim, double> DSvolDM( const FTensors&, const FTensors&, const FTensors&, const double );                       //!< Derivative of the volumetric part of S (Second Piola-Kirchoff stress) wrt to M - see notes
  
  ttl::Tensor<4, dim, double> DPsidDM( const FTensors&, const FTensors&, const FTensors&, const FTensors&, const double );      //!< Derivative of Psid wrt to M - see notes

  FTensors DgammadotvDM(const ttl::Tensor<4, dim, double>&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const double );    //!< Derivative of gammadotv wrt to M - see notes
  FTensors DgammadotvDM(const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const double );    //!< Derivative of gammadotv wrt to M - see notes

  ttl::Tensor<4, dim, double> DSDF(ttl::Tensor<4, dim, double>, const double, const bool);                                       //!< This function evaluates the fourth order tensor DS/DF (total derivative) with dmdf passed as an argument
 
public:
  
  // Constructors
  KMS_IJSS2017_Implicit_BE() : KMS_IJSS2017_Integration_Algorithms<dim> () {}        //!< Default Constructor
  KMS_IJSS2017_Implicit_BE(KMS_IJSS2017_Parameters* , KMS_IJSS2017_IO*, TimeIntegrationDataManager*);              //!< Constructor with material parameters and IO
  
  // Destructors
  virtual ~KMS_IJSS2017_Implicit_BE(){}
  
  // Integrator methods: diagonal terms of the NR matrix
  ttl::Tensor<4, dim, double> DRMDM( const double, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const double );     //!< Diagonal contribution of the residual wrt to M
  ttl::Tensor<4, dim, double> VerboseDRMDM( const double, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const double );     //!< Diagonal contribution of the residual wrt to M
  double DRpcDpc( const double, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const double );                        //!< Diagonal contribution of the residual wrt to pc

  // Integrator methods: right hand side
  ttl::Tensor<2, dim, double> RM( const double, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const double );     //!< Right hand side contribution of the residual for M
  double Rpc( const double, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const FTensors&, const double );                         //!< Right hand side contribution of the residual for pc

  // IO
  void AsAString( std::string&, bool );
  
  // Test methods
  void TestSuite( std::string&, const FTensors&, const FTensors& );                                      //!< This method tests some of the features of the KMS_IJSS2017_Implicit_BE integrator

};



template <int dim>
class KMS_IJSS2017_Implicit_BE_Staggered   : public KMS_IJSS2017_Implicit_BE<dim>
/*!
 \brief This class contains data and methods for the KMS-IJSS2017 model and the algorithms for the implicit integrator
 */
{
  
  typedef ttl::Tensor<2, dim, double> FTensors;

private:
  
  double STAGGEREDTOL=1E-12;
  unsigned MAX_N_OF_STAGGERED_ITERATIONS=100;
  
  void ttlsolveExceptionHandling( const ttl::Tensor<4, dim, double>&, const FTensors& );
  void ttlinverseExceptionHandling( const ttl::Tensor<4, dim, double>& );
  
public:
  
  // Constructors
  KMS_IJSS2017_Implicit_BE_Staggered() : KMS_IJSS2017_Implicit_BE<dim> () {}           //!< Default Constructor
  KMS_IJSS2017_Implicit_BE_Staggered(KMS_IJSS2017_Parameters* , KMS_IJSS2017_IO*, TimeIntegrationDataManager*);      //!< Constructor with material parameters and IO
  
  // Destructors
  virtual ~KMS_IJSS2017_Implicit_BE_Staggered(){}
  
  // IO
  void AsAString( std::string&, bool );
  
  // Integrator methods
  unsigned StepUpdate( const FTensors&, const double, bool );      //!< Updates stresses and internal variables from a given FnP1 and dt
  unsigned VerboseStepUpdate( const FTensors&, const double, bool );      //!< Updates stresses and internal variables from a given FnP1 and dt
  
  unsigned FindpcFromJp( const double, double&, double, bool = false);               //!< Estimate pc at iteration r of the staggered algorithm from Jp at the same step

  // Integrator methods: DMDF
  ttl::Tensor<4, dim, double> IterativeDMDF(ttl::Tensor<4, dim, double>, const double, const bool);                              //!< This function evaluates iteratively the fourth order tensor DM/DF
  ttl::Tensor<4, dim, double> DMDF(const double, const bool);                                                                    //!< This function evaluates the fourth order tensor DM/DF
  void DMDFandDSDF(ttl::Tensor<4, dim, double> &, ttl::Tensor<4, dim, double> &, const double, const bool);                      //!< This function evaluates the fourth order tensors DM/DF and DS/DF (total derivatives)
  

};


template <int dim>
class KMS_IJSS2017_Implicit_BE_Monolithic   : public KMS_IJSS2017_Implicit_BE<dim>
/*!
 \brief This class contains data and methods for the KMS-IJSS2017 model and the algorithms for the implicit integrator
 */
{
  
  typedef ttl::Tensor<2, dim, double> FTensors;
  
public:
  
  // Constructors
  KMS_IJSS2017_Implicit_BE_Monolithic() : KMS_IJSS2017_Implicit_BE<dim> () {}        //!< Default Constructor
  KMS_IJSS2017_Implicit_BE_Monolithic(KMS_IJSS2017_Parameters* , KMS_IJSS2017_IO*);              //!< Constructor with material parameters and IO
  
  // Destructors
  virtual ~KMS_IJSS2017_Implicit_BE_Monolithic(){}
  
  // Integrator methods
  void dRMdpc();     //!< Off-diagonal contribution of the residual wrt to M
  void dRpcdM();     //!< Off-diagonal contribution of the residual wrt to pc
  
  // IO
  void AsAString( std::string&, bool );
};




#include "MechanicalModels_Templates.hpp"
#include "KMS-IJSS2017-templates.hpp"
#include "KMS-IJSS2017-integrator-templates.hpp"
#include "KMS-IJSS2017-explicit-integrator-templates.hpp"
#include "KMS-IJSS2017-implicit-integrator-templates.hpp"
#include "KMS-IJSS2017-implicit-staggered.hpp"


#endif /* KMS_IJSS2017_h */
