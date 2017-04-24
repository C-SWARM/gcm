//
//  TimeIntegration_Manager.cpp
//  ttl-learning
//
//  Created by Alberto Salvadori on 3/17/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#include "TimeIntegration_Manager.h"


// Constructors

TimeIntegrationDataManager::TimeIntegrationDataManager(const double tt0, const double ttf, const double tdt, const double tcurrtime, const size_t tstep)
//! Default constructor
{
  SetTimeFrame( tt0, ttf, tdt, tcurrtime, tstep);
}


// Methods


int TimeIntegrationDataManager::NewStep()
//! Increments the time step. Controls that the time does not exceed the final time,
//! in which case a negative number is returned
{
  timestep++;
  currenttime += deltat;
  
  if ( currenttime > tf )
    return -1;
  
  return 1;
}


void TimeIntegrationDataManager::SetTimeFrame(const double tt0, const double ttf, const double tdt, const double tcurrtime, const size_t tstep)
//! Initialize the integrator with t0, tf, deltat
{
  t0 = tt0;
  tf = ttf;
  deltat = tdt;
  currenttime = tcurrtime;
  timestep = tstep;
}

