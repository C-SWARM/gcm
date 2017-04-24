//
//  MechanicalModels.h
//  ttl-learning
//
//  Created by Alberto Salvadori on 2/1/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//


/*!
 
 \file MechanicalModels.h
 \brief Header file for abstract class of material models
 
 */


#ifndef MechanicalModels_h
#define MechanicalModels_h


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ttl/ttl.h>


class MechanicalModels_IO
/*!
 \brief This class contains and handles the input and output for the KMS-IJSS2017 model
 */
{
protected:
  std::ofstream logfile;    //!< log file for KMS_IJSS2017
  std::ofstream sigmafile;  //!< output of the Cauchy stress for KMS_IJSS2017
  
  size_t PrintEveryNSteps;
  
  
public:
  // Constructors & destructors
  MechanicalModels_IO(){}
  MechanicalModels_IO( const std::string&, const std::string&, unsigned = 1 );
  virtual ~MechanicalModels_IO();
  
  // Accessing IO
  std::ofstream& log(){ return logfile; }
  std::ofstream& sigma(){ return sigmafile; }
  
  bool StepToPrint( size_t );
  
};




template <int dim>
class MechanicalModels_Test_Functions
/*!
 \brief This class contains and handles the test functions for a generic material model
 */
{
  
public:
  
  typedef ttl::Tensor<2, dim, double> FTensors;
  
  // Integration Test Methods
  static void IsotropicCompaction(const FTensors&, double, FTensors&);              //!< Isotropic compaction test function
  static void ReversedIsotropicCompaction(const FTensors&, double, FTensors&);      //!< Isotropic compaction test function, reversed in time
  static void CompactionAndShear1(const FTensors&, double, FTensors&);              //!< Isotropic compaction up to 40s, then simple shear
  static void ReversedCompactionAndShear1(const FTensors&, double, FTensors&);      //!< Isotropic compaction test function, reversed in time
  static void SmallStrainCompactionAndShear2(const FTensors&, double, FTensors&);   //!< Isotropic compaction up to 40s, then simple shear
  static void SmallStrainShear(const FTensors&, double, FTensors&);                 //!< Simple shear
  
};




template <int dim>
class MechanicalModels
/*!
 \brief This is an abstract class for a generic material model
 */
{
  

public:
  
  MechanicalModels(){}
  
  // IO
  void AsAString( std::string&, bool ){}

  
};





#endif /* MechanicalModels_h */
