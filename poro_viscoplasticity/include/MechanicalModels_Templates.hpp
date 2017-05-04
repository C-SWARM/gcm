//
//  MechanicalModels_Templates.hpp
//  ttl-learning
//
//  Created by Alberto Salvadori on 2/2/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#ifndef MechanicalModels_Templates_h
#define MechanicalModels_Templates_h


#include <math.h>
#include <random>

#include "MechanicalModels.h"
#include "ttl-tools.h"




// class MechanicalModels_Test_Functions
// *********************************



// Histories for F(t)
// ------------------
template <int dim>
void MechanicalModels_Test_Functions<dim>::IsotropicCompaction(const FTensors& F0, double t, FTensors& F)
//! This function implements an isotropic, triaxial compaction.
//! Off-diagonal terms are perturbations close to zero.
{
  double strainrate = -0.005, perturbator = 1E-13;
  
  using namespace ttlindexes;
  
  // FTensors Id2 = ttl::identity(i,j);
 
  std::mt19937 gen;
  std::uniform_real_distribution<double> dis(0.0,1.0);
  
  FTensors PerturbedId2 = {  1, dis(gen) * perturbator, dis(gen) * perturbator,
                             dis(gen) * perturbator, 1, dis(gen) * perturbator,
                             dis(gen) * perturbator, dis(gen) * perturbator, 1 };
  
  
  F(i,j) = F0(i,j) + strainrate * t * PerturbedId2(i,j);
}



template <int dim>
void MechanicalModels_Test_Functions<dim>::ReversedIsotropicCompaction(const FTensors& F0, double t, FTensors& F)
{
  
  using namespace ttlindexes;

  double strainrate = -0.005;
  
  FTensors Id2 = ttl::identity(i,j);
  
  F(i,j) = F0(i,j) + strainrate * t * sin( t * M_PI/66.0 ) * Id2(i,j);
}


template <int dim>
void MechanicalModels_Test_Functions<dim>::ReversedCompactionAndShear1(const FTensors& F0, double t, FTensors& F)
//! This function implements an isotropic, triaxial compaction, with sign reversal and shear.
//! Off-diagonal terms are perturbations close to zero.
{
  
  using namespace ttlindexes;
  
    // at t=42.6 pc achieves its maximum in this experiment
  
  double strainrate = -0.005, perturbator = 1E-91;

  std::mt19937 gen;
  std::uniform_real_distribution<double> dis(0.0,1.0);
  
  FTensors PerturbedId2 = {  1, dis(gen) * perturbator, dis(gen) * perturbator,
    dis(gen) * perturbator, 1, dis(gen) * perturbator,
    dis(gen) * perturbator, dis(gen) * perturbator, 1 };
  
  FTensors PerturbedSh2 = {  dis(gen) * perturbator, 1, dis(gen) * perturbator,
    dis(gen) * perturbator, dis(gen) * perturbator, dis(gen) * perturbator,
    dis(gen) * perturbator, dis(gen) * perturbator, dis(gen) * perturbator };

  double lambdac = (t > 46) ?strainrate * 46 * sin( 46 * M_PI/66.0 )  :  strainrate * t * sin( t * M_PI/66.0 )  ;
  double lambdas = (t > 42.6) ? ( 42.6 - t ) * strainrate : 0;
  
  // FTensors Id2 = ttl::identity(i,j);
  // FTensors Sh2 = {0,1,0,0,0,0,0,0,0};
  
  F(i,j) = F0(i,j) + lambdac * PerturbedId2(i,j) + lambdas * PerturbedSh2(i,j);
}



template <int dim>
void MechanicalModels_Test_Functions<dim>::CompactionAndShear1(const FTensors& F0, double t, FTensors& F)
//! History for F, to be used in large strain analyses
{
  
  using namespace ttlindexes;
  
  double lambdac = (t > 40) ? 0.2 : t / 200.00;
  double lambdas = (t > 40) ? ( t - 40 )/ 200.00 : 0;
  
  F = { 1-lambdac, lambdas, 0, 0, 1-lambdac, 0, 0, 0, 1-lambdac  };
}




template <int dim>
void MechanicalModels_Test_Functions<dim>::SmallStrainCompactionAndShear2(const FTensors& F0, double t, FTensors& F)
//! History for displacements, to be used in small strain analyses
{
  
  using namespace ttlindexes;
  
  double lambdac = (t > 0.4) ? 0.002 : t / 200.00;
  double lambdas = (t > 0.4) ? ( t - 0.4 )/ 200.00 : 0;
  
  F = { -lambdac, lambdas, 0, lambdas, -lambdac, 0, 0, 0, -lambdac  };
}



template <int dim>
void MechanicalModels_Test_Functions<dim>::SmallStrainShear(const FTensors& F0, double t, FTensors& F)
//! History for displacements, to be used in small strain analyses
{
  
  using namespace ttlindexes;
  
  double lambdac = 0;
  double lambdas = 0.5 * t;
  
  F = { -lambdac, lambdas, 0, lambdas, -lambdac, 0, 0, 0, -lambdac  };
}














#endif /* MechanicalModels_Templates_h */
