//
//  TimeIntegration_Manager.hpp
//  ttl-learning
//
//  Created by Alberto Salvadori on 3/17/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#ifndef TimeIntegration_Manager_h
#define TimeIntegration_Manager_h

#include <stdio.h>



class TimeIntegrationDataManager
//! Data for time integration.
//! It is shared by all finite-difference time integration schemes,
//! which are assumed to share the same time integration data.
//! In particular, all time steps are the same.
{
protected:
  
  double t0, tf;                                    //!< Initial time, final time
  double deltat;                                    //!< Time increment
  double currenttime;                               //!< Current time
  size_t timestep;                                  //!< Time step
  
  // Integrator methods
  void SetTimeFrame(const double,const double, const double, const double, const size_t);    //!< Initialize the integrator with t0, tf, deltat, currenttime, timestep
  
public:
  
  // Constructors
  TimeIntegrationDataManager()  { t0=tf=deltat=currenttime=0; timestep = 0; };    //!< Default Constructor
  TimeIntegrationDataManager(const double,const double, const double, const double, const size_t);    //!< Default Constructor
  
  // Access data
  double InitialTime(){ return t0; };
  double FinalTime(){ return tf; };
  double Deltat(){ return deltat; };
  double CurrentTime(){ return currenttime; };
  size_t TimeStep(){ return timestep; };
  
  // Increase the time step by one
  int NewStep();
};




#endif /* TimeIntegration_Manager_h */
