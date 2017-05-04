//
//  MechanicalModels_integrator.h
//  ttl-learning
//
//  Created by Alberto Salvadori on 1/25/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#ifndef MechanicalModels_Integrator_h
#define MechanicalModels_Integrator_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ttl/ttl.h>

#include "ttl-tools.h"
#include "TimeIntegration_Manager.h"


namespace DualityPairs
{
  enum Pairs { FirstPiola, SecondPiolaKirchoff, Kirchoff };
}



template <int dim>
class MechanicalModels_FiniteDifference_Integrator
/*!
 \brief This abstract class contains data and methods for the integration of a generic model  */
{
  
  typedef ttl::Tensor<2, dim, double> FTensors;
  typedef ttl::Tensor<1, dim, double> Vectors;
  
protected:
  
  // Data
  double INTEGRATOR_TOL=1E-12;                      //!< Tolerance in this class
  TimeIntegrationDataManager* TimeIntegrationData;  //!< Data for time integration. It is a pointer because all integrators are assumed to share the same time integration data.
  
  // Integrator methods
  void SetTimeFrame(const TimeIntegrationDataManager*);    //!< Initialize the integrator with t0, tf, deltat, currenttime, timestep
  void CleanFTensors( FTensors& );   //!< Sets very small off-diagonal terms to zero
  
public:
  
  // Constructors
  MechanicalModels_FiniteDifference_Integrator()  { }    //!< Default Constructor
  MechanicalModels_FiniteDifference_Integrator( const TimeIntegrationDataManager* tmd )  { TimeIntegrationData = tmd; }    //!< Default Constructor
  
  // Destructors
  virtual ~MechanicalModels_FiniteDifference_Integrator() { }    //!< Destructor

  // Reading methods
  virtual FTensors get_epse(){using namespace ttlindexes; return ttl::zero(i,j);}
  virtual FTensors get_epsp(){using namespace ttlindexes; return ttl::zero(i,j);}
  virtual FTensors get_Cauchy_stress(){using namespace ttlindexes; return ttl::zero(i,j);}
  virtual FTensors get_dev_Cauchy_stress(){using namespace ttlindexes; return ttl::zero(i,j);}
  virtual ttl::Tensor<2, 3, double> get_dev3D_Cauchy_stress(){using namespace ttlindexes; return ttl::zero(i,j);}
  virtual double get_Cauchy_pressure(){return 0;} //!< Cauchy pressure
  virtual double get_Mises_stress(){return 0;} //!< Mises stress
  
  // Integrator methods
  virtual unsigned StepUpdate( const FTensors&, const double, bool );     //!< Updates stresses and internal variables from a given FnP1 and dt

  // ConsistentTangentStiffness methods:
  // These methods return the Consisten Tangent Stiffness matrix
  // and change with the duality pair that has been implemented ( Piola, Second Piola-Kirchoff, Kirchoff, Cauchy)
  // the arguments that are passed are:
  // - the deformation measure at previous iteration at the Gauss point of interest
  // - double deltat
  // - int DualityPair
  // - bool Verbose
  virtual ttl::Tensor<4, dim, double> ConsistentTangentStiffness( const FTensors&, const double = 0, const bool = false );
  virtual ttl::Tensor<4, dim, double> ConsistentTangentStiffness( const double = 0, const bool = false );
  virtual ttl::Tensor<4, dim, double> ConsistentTangentStiffness( const bool = false );
  virtual ttl::Tensor<4, dim, double> ConsistentTangentStiffness( const FTensors&, const double = 0, const int = 0, const bool = false );
  virtual ttl::Tensor<4, dim, double> ConsistentTangentStiffness( const double = 0, const int = 0, const bool = false );

  virtual ttl::Tensor<4, 3, double> ConsistentTangentStiffness_for_dev( const double = 0, const bool = false );
  virtual double ConsistentTangentStiffness_for_pressure( const double = 0, const bool = false );

  
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
  virtual FTensors GateauxDifferential( const ttl::Tensor<4, dim, double>&, const FTensors&, const Vectors&, double, int = 0, bool = false );
  virtual FTensors GateauxDifferential( const ttl::Tensor<4, dim, double>&, const FTensors&, double, int = 0, bool = false );
  virtual FTensors GateauxDifferential( const FTensors&, const FTensors&, double = 0, int = 0, bool = false );
  virtual FTensors GateauxDifferential( const FTensors&, const Vectors&, double, int = 0, bool = false );
  virtual FTensors GateauxDifferential( const FTensors&, double = 0, int = 0, bool = false );
};


/*
template <int dim>
void MechanicalModels_FiniteDifference_Integrator<dim>::SetTimeFrame(const double tt0, const double ttf, const double tdt, const double tcurrtime, const size_t tstep)
//! Initialize the integrator with t0, tf, deltat
{
  t0 = tt0;
  tf = ttf;
  Deltat = tdt;
  CurrentTime = tcurrtime;
  timestep = tstep;
}
*/


template <int dim>
void MechanicalModels_FiniteDifference_Integrator<dim>::SetTimeFrame(const TimeIntegrationDataManager* tdm)
//!< Initialize the integrator with the TimeIntegrationDataManager
{
  TimeIntegrationData = tdm;
}


template <int dim>
void MechanicalModels_FiniteDifference_Integrator<dim>::CleanFTensors( FTensors& WhateverF )
//! Sets very small off-diagonal terms to zero
//! The need of this method stands in the fact that very small off-diagonal terms (of the order E-90)
//! cause issues in the Mechanical Model Integration. Rather, if they are set to zero, no issues whatsoever
//! arise. The check can be done on F since the diagonal term should be of magnitude 1.
{
  double FrNorm = FrobeniusNorm( WhateverF );
  double ZERO = FrNorm * 1E-16;
  
  for ( unsigned i=0; i<dim; i++ )
    for ( unsigned j=0; j<dim; j++ )
      if ( fabs( WhateverF[i][j] ) < ZERO )
        WhateverF[i][j] = 0.0 ;
  
}


template <int dim>
unsigned MechanicalModels_FiniteDifference_Integrator<dim>::StepUpdate( const FTensors&, const double, bool )
//! Updates stresses and internal variables from a given FnP1 and dt
{ return 1; }


// ConsistentTangentStiffness methods:
// These methods return the Consistent Tangent Stiffness matrix
// and change with the duality pair that has been implemented ( Piola, Second Piola-Kirchoff, Kirchoff, Cauchy)
// the arguments that are passed are:
// - the deformation measure at previous iteration at the Gauss point of interest
// - double deltat
// - int DualityPair
// - bool Verbose
template <int dim>
ttl::Tensor<4, dim, double> MechanicalModels_FiniteDifference_Integrator<dim>::
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
ttl::Tensor<4, dim, double> MechanicalModels_FiniteDifference_Integrator<dim>::ConsistentTangentStiffness(
                                                                                                          const double dt,
                                                                                                          const int DualityPair,
                                                                                                          const bool Verbose
                                                                                                          )
//! Updates the tangent stiffness matrix from a given FnP1 and dt
{
  return ConsistentTangentStiffness( Verbose );
}

template <int dim>
ttl::Tensor<4, dim, double> MechanicalModels_FiniteDifference_Integrator<dim>::
ConsistentTangentStiffness(
                                                               const FTensors& FnP1,
                                                               const double dt,
                                                               const bool Verbose
                                                               )
//! Updates the tangent stiffness matrix from a given FnP1 and dt
//! This form is tipycally used in small strains, where no issues on the
//! duality pair comes into play
{
  return ConsistentTangentStiffness( Verbose );
}

template <int dim>
ttl::Tensor<4, dim, double> MechanicalModels_FiniteDifference_Integrator<dim>::
ConsistentTangentStiffness(
                           const double dt,
                           const bool Verbose
                           )
//! Updates the tangent stiffness matrix from a given FnP1 and dt
//! This form is tipycally used in small strains, where no issues on the
//! duality pair comes into play
{
  return ConsistentTangentStiffness( Verbose );
}


template <int dim>
ttl::Tensor<4, dim, double> MechanicalModels_FiniteDifference_Integrator<dim>::
ConsistentTangentStiffness( const bool )
//! Updates the tangent stiffness matrix from a given FnP1 and dt
//! This form is tipycally used in small strains, where no issues on the
//! duality pair comes into play
{
  using namespace ttlindexes;
  
  ttl::Tensor<4, dim, double> Id4 = ttl::identity<dim>( i,j,k,l );
  return Id4;
}


template <int dim>
ttl::Tensor<4, 3, double> MechanicalModels_FiniteDifference_Integrator<dim>::
ConsistentTangentStiffness_for_dev(
                                   const double dt,
                                   const bool Verbose
                                   )
//! Updates the tangent stiffness matrix for the deviatoric, isochoric part
//! This form is tipycally used in small strains, where no issues on the
//! duality pair comes into play
//! ConsistentTangentStiffness_for_dev is always in dimension 3
//1 See ttl-tools.h
{
  using namespace ttlindexes;
  ttl::Tensor<4, 3, double> Id4 = ttl::identity<3>( i,j,k,l );
  return Id4;
}


template <int dim>
double MechanicalModels_FiniteDifference_Integrator<dim>::
ConsistentTangentStiffness_for_pressure(
                                        const double dt,
                                        const bool Verbose
                                        )
//! Updates the tangent stiffness matrix for the volumetric part
//! This form is tipycally used in small strains, where no issues on the
//! duality pair comes into play
{
  return -1;
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
ttl::Tensor<2, dim, double> MechanicalModels_FiniteDifference_Integrator<dim>::
      GateauxDifferential(
                          const ttl::Tensor<4, dim, double>& TangStiff,
                          const FTensors& shfgrad,
                          const Vectors& sol,
                          double dt,
                          int DualityPair,
                          bool Verbose
                          )
{
  return GateauxDifferential( shfgrad, sol, dt, DualityPair, Verbose );
}

template <int dim>
ttl::Tensor<2, dim, double> MechanicalModels_FiniteDifference_Integrator<dim>::
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
ttl::Tensor<2, dim, double> MechanicalModels_FiniteDifference_Integrator<dim>::
GateauxDifferential(
                    const FTensors& shfgrad,
                    const Vectors& sol,
                    double dt,
                    int DualityPair,
                    bool Verbose
                    )
{
  return GateauxDifferential( shfgrad, dt, DualityPair, Verbose  );
}

template <int dim>
ttl::Tensor<2, dim, double> MechanicalModels_FiniteDifference_Integrator<dim>::
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
ttl::Tensor<2, dim, double> MechanicalModels_FiniteDifference_Integrator<dim>::
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





#endif /* MechanicalModels_Integrator_h */
