//
//  FeFp_integrator_templates.hpp
//  ttl-learning
//
//  Created by alberto salvadori on 3/1/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#ifndef FeFp_integrator_templates_h
#define FeFp_integrator_templates_h




// Protected Methods
// -----------------

template <int dim>
void FeFp_Integrator<dim>::StepClose( )
//void KMS_IJSS2017_Integrator<dim>::StepClose( )
//! Assigns data at step n from data at step n+1
{
  
  FiniteStrainModels_Integrator<dim>::StepClose( );
  
  eFn = eFnp1 ;              // Elastic part deformation gradient F
  pFn = pFnp1;               // Plastic part deformation gradient F
  
}


// Public Methods
// --------------



template <int dim>
unsigned FeFp_Integrator<dim>::StepUpdate( const FTensors&, const double, bool )
//! Updates stresses and internal variables from a given FnP1 and dt
{ return 1; }






template <int dim>
ttl::Tensor<4, dim, double> FeFp_Integrator<dim>::DMDF(
                                                       const double dt,
                                                       const bool Verbose
                                                       )
//! This function evaluates the fourth order tensor DM/DF and stores it into dmdf
//! It depends upon every constitutive model, therefore the function is declared as
//! virtual
{
  using namespace ttlindexes;
  return ttl::identity<dim>( i,j,k,l );
}

template <int dim>
ttl::Tensor<4, dim, double> FeFp_Integrator<dim>::DSDF(
                                                       const double dt,
                                                       const bool Verbose
                                                       )
//! This function evaluates the fourth order tensor DS/DF
//! It depends upon every constitutive model, therefore the function is declared as
//! virtual
{
  using namespace ttlindexes;
  return ttl::identity<dim>( i,j,k,l );
}

template <int dim>
ttl::Tensor<4, dim, double> FeFp_Integrator<dim>::DSDF(
                                                       ttl::Tensor<4, dim, double> localdmdf,
                                                       const double dt,
                                                       const bool Verbose
                                                       )
//! This function evaluates the fourth order tensor DS/DF using dmdf as an argument that is passed
//! It depends upon every constitutive model, therefore the function is declared as
//! virtual
{
  return DSDF( dt, Verbose );
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
ttl::Tensor<4, dim, double> FeFp_Integrator<dim>::
ConsistentTangentStiffness(
                          const FTensors& FnP1,
                          const double dt,
                          const int DualityPair,
                          const bool Verbose
                          )
//! Updates the tangent stiffness matrix from a given FnP1 and dt
//! As first, since Fnp1 is given, the Step is updated. Then, the constistent tangent stiffness is invoked.
//! Since however it is expected that the Step will be updated elsewhere, this function may actually
//! never be called.
{
  // Update step
  StepUpdate( FnP1, dt, Verbose );
  
  // Constant Tangent Stiffness evaluation
  return ConsistentTangentStiffness( dt, DualityPair, Verbose );
  
}

template <int dim>
ttl::Tensor<4, dim, double> FeFp_Integrator<dim>::
ConsistentTangentStiffness(
                          const double dt,
                          const int DualityPair,
                          const bool Verbose
                          )
//! Updates the tangent stiffness matrix from a given FnP1 and dt
//! This method expects that the Step has already been updated elsewhere, since it uses quantities at step n+1
//! that MUST have been updated before
//!
//! Based on this definition of the ConsistentTangentStiffness, every FeFp - based constitutive model can be integated
//! provided that the functions StepUpdate, DMDF, and DSDF are given
{
  
  ttl::Tensor<4, dim, double> CTS = {};

  switch( DualityPair )
  {
      
    case DualityPairs::FirstPiola :
    {
      using namespace ttlindexes;
      static constexpr ttl::Index<'b'> b;
      static constexpr ttl::Index<'c'> c;
      static constexpr ttl::Index<'d'> d;
      static constexpr ttl::Index<'g'> g;
      
      // calculate DMDF and store in the tensor dmdf
      ttl::Tensor<4, dim, double> dmdf = DMDF( dt, Verbose );
      
      // calculate DSDF and store in the tensor dsdf
      ttl::Tensor<4, dim, double> dsdf =  DSDF( dmdf, dt, Verbose );
      
      // Invert tensor pFnp1
      FTensors InvpFnp1 = ttl::inverse( this->pFnp1 );
      FTensors InvpFn = ttl::inverse( this->pFn );
      
      
      CTS( i,j,k,l ) =
      ttl::identity<dim>( i,b,k,l ) * InvpFnp1( b,c ) * this->Snp1( c,d ) * InvpFnp1( j,d ) +
      this->Fnp1( i,b ) * InvpFn( b,g ) * dmdf( g,c,k,l ) * this->Snp1( c,d ) * InvpFnp1( j,d ) +
      this->Fnp1( i,b ) * InvpFnp1( b,c ) * dsdf( c,d,k,l ) * InvpFnp1( j,d ) +
      this->Fnp1( i,b ) * InvpFnp1( b,c ) * this->Snp1( c,d ) * InvpFn( j,g ) * dmdf( g,d,k,l ) ;
      
      break;
    }
      
      std::cout << " FATAL ERROR: FeFp_Integrator<dim>::ConsistentTangentStiffness( const double dt, const int DualityPair, const bool Verbose) - Unknown duality pair.";
      std::cout << " Code terminates. \n";
      std::exit(EXIT_FAILURE);
  
  }
  
  return CTS;
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
ttl::Tensor<2, dim, double> FeFp_Integrator<dim>::
GateauxDifferential(
                    const ttl::Tensor<4, dim, double>& TangStiff,
                    const FTensors& shfgrad,
                    double dt,
                    int DualityPair,
                    bool Verbose
                    )
//! Generic GateauxDifferential, merely defined as the product of the tangent stiffness and
//! the shape functions gradient tensor
{
  // pre-defined ttlindexes from i to l
  using namespace ttlindexes;
  
  // ttl::Tensor<4, dim, double> stiffness;
  ttl::Tensor<2, dim, double> GD;
  
  GD( i,j ) = TangStiff( i,j,k,l ) * shfgrad( k,l );
  
  return GD;
}

template <int dim>
ttl::Tensor<2, dim, double> FeFp_Integrator<dim>::
GateauxDifferential(
                    const FTensors& shfgrad,
                    const FTensors& FnP1,
                    double dt,
                    int DualityPair,
                    bool Verbose
                    )
//! GateauxDifferential which constructs the consistent tangent stiffness
//! at the give Gauss point before multiplication.
//! Since this function is called with the tensor FnP1, it is assumed
//! that the Step has not been yet updated and the ConsistentTangentStiffness method
//! will take care of that.
{
  // pre-defined ttlindexes from i to l
  using namespace ttlindexes;
  
  // estimation of the tangent stiffness, including Step Update
  ttl::Tensor<4, dim, double> TangStiff = ConsistentTangentStiffness( FnP1, dt, DualityPair, Verbose );
  
  // ttl::Tensor<4, dim, double> stiffness;
  ttl::Tensor<2, dim, double> GD;
  GD( i,j ) = TangStiff( i,j,k,l ) * shfgrad( k,l );
  
  return GD;
}

template <int dim>
ttl::Tensor<2, dim, double> FeFp_Integrator<dim>::
GateauxDifferential(
                    const FTensors& shfgrad,
                    double dt,
                    int DualityPair,
                    bool Verbose
                    )
//! GateauxDifferential which constructs the consistent tangent stiffness
//! at the give Gauss point before multiplication.
//! Since this function is called without the tensor FnP1, it is assumed
//! that the Step has been updated elsewhere and the ConsistentTangentStiffness method
//! will not take care of that.
{
  // pre-defined ttlindexes from i to l
  using namespace ttlindexes;
  
  // estimation of the tangent stiffness, including Step Update
  ttl::Tensor<4, dim, double> TangStiff = ConsistentTangentStiffness( dt, DualityPair, Verbose );
  
  // ttl::Tensor<4, dim, double> stiffness;
  ttl::Tensor<2, dim, double> GD;
  GD( i,j ) = TangStiff( i,j,k,l ) * shfgrad( k,l );
  
  return GD;
  
}





#endif /* FeFp_integrator_templates_h */
