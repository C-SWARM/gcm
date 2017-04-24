//
//  FeFp_integrator.h
//  ttl-learning
//
//  Created by Alberto Salvadori on 1/25/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#ifndef FeFp_integrator_h
#define FeFp_integrator_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ttl/ttl.h>

#include "FiniteStrainModels_integrator.h"
#include "SmallStrainIntegrator.h"
#include "ttl-tools.h"




// class FeFp_Integrator
// *********************

template <int dim>
class FeFp_Integrator : public FiniteStrainModels_Integrator<dim>
/*!
 \brief This abstract class contains data and methods for the integration of a generic model that
 is based on the FeFp decomposition (or analogous ones, as FeFt)
 */
{
  
  typedef ttl::Tensor<2, dim, double> FTensors;
  typedef ttl::Tensor<1, dim, double> Vectors;
  
protected:
  
  FTensors eFn;                                     //!< Elastic part deformation gradient F at step n. Passed from the external NR algorithm
  FTensors pFn;                                     //!< Plastic part deformation gradient F at step n. Passed from the external NR algorithm
  FTensors eFnp1;                                   //!< Elastic part deformation gradient F at step n+1
  FTensors pFnp1;                                   //!< Plastic part deformation gradient F at step n+1
  FTensors Dpn;                                     //!< Plastic stretching Dp at step n.
  
  // Integrator methods
  void StepClose( );                                //!< Assigns data at step n from data at step n+1
  
public:
  
  // Constructors
  FeFp_Integrator()  { }          //!< Default Constructor
  
  // Integrator methods
  virtual unsigned StepUpdate( const FTensors&, const double, bool );                               //!< Updates stresses and internal variables from a given FnP1 and dt

  virtual ttl::Tensor<4, dim, double> DMDF(const double, const bool);                               //!< This function evaluates the fourth order tensor DM/DF
  virtual ttl::Tensor<4, dim, double> DSDF(const double, const bool);                               //!< This function evaluates the fourth order tensor DS/DF (total derivative)
  virtual ttl::Tensor<4, dim, double> DSDF(ttl::Tensor<4, dim, double>, const double, const bool);  //!< This function evaluates the fourth order tensor DS/DF with dmdf passed as an argument
  

  // ConsistentTangentStiffness methods:
  // These methods return the Consisten Tangent Stiffness matrix
  // and change with the duality pair that has been implemented ( Piola, Second Piola-Kirchoff, Kirchoff, Cauchy)
  // the arguments that are passed are:
  // - the deformation measure at previous iteration at the Gauss point of interest
  // - double deltat
  // - int DualityPair
  // - bool Verbose
  virtual ttl::Tensor<4, dim, double> ConsistentTangentStiffness( const double, const int = 0, const bool Verbose = false );
  virtual ttl::Tensor<4, dim, double> ConsistentTangentStiffness( const FTensors&, const double, const int = 0, const bool Verbose = false );
  
  // Gateaux differential methods:
  // These methods return the "product" of the ConsistentTangentStiffness and the shape function matrix
  // and change with the duality pair that has been implemented ( Piola, Second Piola-Kirchoff, Kirchoff, Cauchy)
  // the arguments that are passed are:
  // - ttl::Tensor<4, dim, double> ConsistentTangentStiffness ( when constant ... )
  // - FTensors shape function matrix
  // - Vector solution at previous iteration at the Gauss point of interest
  // - double deltat
  // - int DualityPair
  // - bool Verbose
  virtual FTensors GateauxDifferential( const ttl::Tensor<4, dim, double>&, const FTensors&,  double = 0, int = 0, bool Verbose = false );
  virtual FTensors GateauxDifferential( const FTensors&, const FTensors&, double = 0, int = 0, bool Verbose = false );
  virtual FTensors GateauxDifferential( const FTensors&, double = 0, int = 0, bool Verbose = false );
  
};

#include "FeFp_integrator_templates.hpp"

#endif /* FeFp_integrator_h */
