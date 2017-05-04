//
//  FiniteStrainModels_integrator.h
//  ttl-learning
//
//  Created by Alberto Salvadori on 1/25/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#ifndef FiniteStrainModels_integrator_h
#define FiniteStrainModels_integrator_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ttl/ttl.h>

#include "MechanicalModels_Integrator.h"
#include "ttl-tools.h"



template <int dim>
class FiniteStrainModels_Integrator : public MechanicalModels_FiniteDifference_Integrator<dim>
/*!
 \brief This abstract class contains data and methods for the integration of a generic model  */
{
  
  typedef ttl::Tensor<2, dim, double> FTensors;
  typedef ttl::Tensor<1, dim, double> Vectors;
  
protected:
  
  // Data
  FTensors Fn;                                      //!< Total deformation gradient F at step n. Passed from the external NR algorithm
  FTensors Fnp1;                                    //!< Total deformation gradient F at step n+1. Passed from the external NR algorithm
  FTensors Sn;                                      //!< Second Piola-Kirchoff stress tensor S at step n.
  FTensors Snp1;                                    //!< Second Piola-Kirchoff stress tensor S at step n+1.
  FTensors KSn;                                     //!< Kirchoff stress tensor S at step n.
  FTensors KSnp1;                                   //!< Kirchoff stress tensor S at step n+1.
  FTensors sigman;                                  //!< Cauchy stress tensor S at step n.
  FTensors sigmanp1;                                //!< Cauchy stress tensor S at step n+1.
  
  // Integrator methods
  void StepClose( );                                                    //!< Assigns data at step n from data at step n+1
  
public:
  
  // Constructors
  FiniteStrainModels_Integrator() : MechanicalModels_FiniteDifference_Integrator<dim>(){}    //!< Default Constructor
  
  // Integrator methods
  virtual unsigned StepUpdate( const FTensors&, const double, bool );                                       //!< Updates stresses and internal variables from a given FnP1 and dt

  // ConsistentTangentStiffness methods:
  // These methods return the Consisten Tangent Stiffness matrix
  // and change with the duality pair that has been implemented ( Piola, Second Piola-Kirchoff, Kirchoff, Cauchy)
  // the arguments that are passed are:
  // - the deformation measure at previous iteration at the Gauss point of interest
  // - double deltat
  // - int DualityPair
  // - bool Verbose
  virtual ttl::Tensor<4, dim, double> ConsistentTangentStiffness( const FTensors&, const double = 0, const int = 0, const bool = false );
  virtual ttl::Tensor<4, dim, double> ConsistentTangentStiffness( const double = 0, const int = 0, const bool = false );
  
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
  virtual FTensors GateauxDifferential( const ttl::Tensor<4, dim, double>&, const FTensors&,  double = 0, int = 0, bool  = false );
  virtual FTensors GateauxDifferential( const FTensors&, const FTensors&, double = 0, int = 0, bool  = false );
  virtual FTensors GateauxDifferential( const FTensors&, double = 0, int = 0, bool  = false );
  
};


template <int dim>
unsigned FiniteStrainModels_Integrator<dim>::StepUpdate( const FTensors&, const double, bool )
//! Updates stresses and internal variables from a given FnP1 and dt
{ return 1; }



template <int dim>
void FiniteStrainModels_Integrator<dim>::StepClose( )
//! Assigns data at step n from data at step n+1
{
  Fn = Fnp1;                 // Total deformation gradient F
  Sn = Snp1;                 // Second Piola-Kirchoff stress tensor S
  KSn = KSnp1;               // Kirchoff stress tensor KS
  sigman = sigmanp1;         // Cauchy stress tensor sigma
}



// ConsistentTangentStiffness methods:
// These methods return the Consisten Tangent Stiffness matrix
// and change with the duality pair that has been implemented ( Piola, Second Piola-Kirchoff, Kirchoff, Cauchy)
// the arguments that are passed are:
// - the deformation measure at previous iteration at the Gauss point of interest
// - double deltat
// - int DualityPair
// - bool Verbose
template <int dim>
ttl::Tensor<4, dim, double> FiniteStrainModels_Integrator<dim>::
ConsistentTangentStiffness(
                          const FTensors& FnP1,
                          const double dt,
                          const int DualityPair,
                          const bool Verbose
                          )
//! Updates the tangent stiffness matrix from a given FnP1 and dt
{
  return ConsistentTangentStiffness( dt, DualityPair, Verbose );
}

template <int dim>
ttl::Tensor<4, dim, double> FiniteStrainModels_Integrator<dim>::ConsistentTangentStiffness( const double, const int, const bool )
//! Updates the tangent stiffness matrix from a given FnP1 and dt
{
  using namespace ttlindexes;
  
  ttl::Tensor<4, dim, double> Id4 = ttl::identity<dim>( i,j,k,l );
  return Id4;
}



// Gateaux differential methods:
// These methods return the "product" of the ConsistentTangentStiffness and the shape function matrix
// and change with the duality pair that has been implemented ( Piola, Second Piola-Kirchoff, Kirchoff, Cauchy)
// the arguments that are passed are:
// - int, DualityPair
// - ttl::Tensor<4, dim, double> ConsistentTangentStiffness ( when constant ... )
// - FTensors shape function matrix
// - Vector solution at previous iteration at the Gauss point of interest

template <int dim>
ttl::Tensor<2, dim, double> FiniteStrainModels_Integrator<dim>::
GateauxDifferential(
                    const ttl::Tensor<4, dim, double>& TangStiff,
                    const FTensors& shfgrad,
                    double dt,
                    int DualityPair,
                    bool Verbose
                    )
{
  return GateauxDifferential( shfgrad, dt, DualityPair, Verbose );
}

template <int dim>
ttl::Tensor<2, dim, double> FiniteStrainModels_Integrator<dim>::
GateauxDifferential(
                    const FTensors& shfgrad,
                    const FTensors& Fnp1,
                    double dt,
                    int DualityPair,
                    bool Verbose
                    )
{
  return GateauxDifferential( shfgrad, dt, DualityPair, Verbose  );
}

template <int dim>
ttl::Tensor<2, dim, double> FiniteStrainModels_Integrator<dim>::
GateauxDifferential(
                    const FTensors& shfgrad,
                    double dt,
                    int DualityPair,
                    bool Verbose
                    )
{
  using namespace ttlindexes;
  
  ttl::Tensor<2, dim, double> Id2 = ttl::identity<dim>( i,j );
  return Id2;
}





#endif /* FiniteStrainModels_integrator_h */
