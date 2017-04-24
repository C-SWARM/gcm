//
//  main.cpp
//  ttl-learning
//
//  Created by alberto salvadori on 12/8/16.
//  Copyright © 2016 alberto salvadori. All rights reserved.
//

// include

#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

#include <ttl/ttl.h>

#include "KMS-IJSS2017.h"
#include "ttl-tools.h"



void KMS_IJSS2017_Explicit_test()
{
  
  // data for KMS_IJSS2017 definition
  const int spacedim = 3;
  std::string outstr = "";
  
  std::string pathstr = "./out/KMS_IJSS2017";
    
  clock_t t;
  
  unsigned PrintEveryNSteps= 100;
  
  {
    std::string logstr = pathstr + ".explicit.smb.log";
    std::string Fstr = pathstr + ".explicit.smb.F.txt";
    std::string pFstr = pathstr + ".explicit.smb.Fp.txt";
    std::string Sstr = pathstr + ".explicit.smb.S.txt";
    std::string KSstr = pathstr + ".explicit.smb.KS.txt";
    std::string sigmastr = pathstr + ".explicit.smb.sigma.txt";
    std::string pcstr = pathstr + ".explicit.smb.pc.txt";
    std::string matstr = pathstr + ".explicit.smb.mat.txt";
    
    // Time integrator
    double Dt = 0.01, InitialTime=0, FinalTime = 13200 * Dt, CurrentTime=0;
    size_t TimeStep = 0;
    TimeIntegrationDataManager TimeIntegrationData( InitialTime, FinalTime, Dt, CurrentTime, TimeStep );

    // Parameters and Model definition
    bool usingSmoothMacauleyBrackets = true;
    KMS_IJSS2017_Parameters matmodel( 1.0, 1.1, 0.15, 0.0005, 0.62, 0.37, 77.22, 13.01, 15, 0.01, 0.2, 5.8, 30, 60, 0.063, 0.008, 2, 1, 290, usingSmoothMacauleyBrackets );
    KMS_IJSS2017_IO IO(logstr, Fstr, pFstr, Sstr, KSstr, sigmastr, pcstr, matstr, PrintEveryNSteps);
    MechanicalModels_Test_Functions<spacedim> KMS_IJSS2017TestFunctions;
    
    KMS_IJSS2017_Explicit_FE<spacedim> ExplicitModelInstance( &matmodel, &IO, &TimeIntegrationData );
    
    // Integrator test
    bool Verbose = false;
    
    t = clock();
    IO.log() << "  ExplicitModelInstance.IntegratorTest for ReversedIsotropicCompaction running. \n";
    ExplicitModelInstance.IntegratorTest( KMS_IJSS2017TestFunctions.ReversedIsotropicCompaction, Verbose );
    t = clock() - t;
    IO.log() << "  ExplicitModelInstance.IntegratorTest for ReversedIsotropicCompaction using Smooth Macauley Brackets terminated in " << ((float)t)/CLOCKS_PER_SEC << "seconds. \n\n";
  }
  /*
  {
    std::string logstr = pathstr + ".explicit.mb.log";
    std::string Fstr = pathstr + ".explicit.mb.F.txt";
    std::string pFstr = pathstr + ".explicit.mb.Fp.txt";
    std::string Sstr = pathstr + ".explicit.mb.S.txt";
    std::string KSstr = pathstr + ".explicit.mb.KS.txt";
    std::string sigmastr = pathstr + ".explicit.mb.sigma.txt";
    std::string pcstr = pathstr + ".explicit.mb.pc.txt";
    std::string matstr = pathstr + ".explicit.mb.mat.txt";

    // Time integrator
    double Dt = 0.1, InitialTime=0, FinalTime = 100 * Dt, CurrentTime=0;
    size_t TimeStep = 0;
    TimeIntegrationDataManager TimeIntegrationData( InitialTime, FinalTime, Dt, CurrentTime, TimeStep );

    // Parameters and Model definition
    bool usingSmoothMacauleyBrackets = false;
    KMS_IJSS2017_Parameters matmodel( 1.0, 1.1, 0.15, 0.0005, 0.62, 0.37, 77.22, 13.01, 15, 0.01, 0.2, 5.8, 30, 60, 0.063, 0.008, 2, 1, 290, usingSmoothMacauleyBrackets) ;
    KMS_IJSS2017_IO IO(logstr, Fstr, pFstr, Sstr, KSstr,sigmastr, pcstr, matstr, PrintEveryNSteps);
    MechanicalModels_Test_Functions<spacedim> KMS_IJSS2017TestFunctions;
    
    KMS_IJSS2017_Explicit_FE<spacedim> ExplicitModelInstance( &matmodel, &IO, &TimeIntegrationData );
    
    // Integrator test
    bool Verbose = false;
    
    t = clock();
    IO.log() << "  ExplicitModelInstance.IntegratorTest for CompactionAndShear1 running. \n";
    ExplicitModelInstance.IntegratorTest( KMS_IJSS2017TestFunctions.CompactionAndShear1, Verbose );
    t = clock() - t;
    IO.log() << "  ExplicitModelInstance.IntegratorTest for CompactionAndShear1 using standard Macauley Brackets terminated in " << ((float)t)/CLOCKS_PER_SEC << "seconds. \n";
  }
  */
}


void KMS_IJSS2017_Implicit_test()
{

  // data for KMS_IJSS2017 definition
  const int spacedim = 3;
  std::string outstr = "";
  
  std::string pathstr = "out/KMS_IJSS2017";
  
  clock_t t;
  
  unsigned PrintEveryNSteps= 100;

  std::string logstr = pathstr + ".implicit.smb.log";
  std::string Fstr = pathstr + ".implicit.smb.F.txt";
  std::string pFstr = pathstr + ".implicit.smb.Fp.txt";
  std::string Sstr = pathstr + ".implicit.smb.S.txt";
  std::string KSstr = pathstr + ".implicit.smb.KS.txt";
  std::string sigmastr = pathstr + ".implicit.smb.sigma.txt";
  std::string pcstr = pathstr + ".implicit.smb.pc.txt";
  std::string matstr = pathstr + ".implicit.smb.mat.txt";
  
  // Time integrator
  double Dt = 0.01, InitialTime=0, FinalTime = 13200 * Dt, CurrentTime=0;
  size_t TimeStep = 0;
  TimeIntegrationDataManager TimeIntegrationData( InitialTime, FinalTime, Dt, CurrentTime, TimeStep );

  // Parameters and Model definition
  bool usingSmoothMacauleyBrackets = true;
  KMS_IJSS2017_Parameters matmodel( 1.0, 1.1, 0.15, 0.0005, 0.62, 0.37, 77.22, 13.01, 15, 0.01, 0.2, 5.8, 30, 60, 0.063, 0.008, 2, 1, 290, usingSmoothMacauleyBrackets ) ;
  KMS_IJSS2017_IO IO(logstr, Fstr, pFstr, Sstr, KSstr, sigmastr, pcstr, matstr, PrintEveryNSteps);
  
  // Loading case
  MechanicalModels_Test_Functions<spacedim> KMS_IJSS2017TestFunctions;
  
  // Instance for testing
  KMS_IJSS2017_Implicit_BE_Staggered<spacedim> ImplicitModelInstance( &matmodel, &IO, &TimeIntegrationData );
  
  // Integrator test
  bool Verbose = false;
  
  t = clock();
  IO.log()  << "\n  ImplicitModelInstance.IntegratorTest for ReversedCompactionAndShear1 running. \n";
  ImplicitModelInstance.IntegratorTest( KMS_IJSS2017TestFunctions.ReversedCompactionAndShear1, Verbose );
  t = clock() - t;
  IO.log()  << "\n  ImplicitModelInstance.IntegratorTest for ReversedCompactionAndShear1 using smooth Macauley Brackets terminated in " << ((float)t)/CLOCKS_PER_SEC << "seconds. \n";
  
}





void KMS_IJSS2017_Implicit_ConsistentTangent_test()
{
  
  // data for KMS_IJSS2017 definition
  const int spacedim = 3;
  std::string outstr = "";
  
  std::string pathstr = "out/KMS_IJSS2017";
  
  clock_t t;
  
  unsigned PrintEveryNSteps= 1;
  
  std::string logstr = pathstr + ".implicit.cst.smb.log";
  std::string Fstr = pathstr + ".implicit.cst.smb.F.txt";
  std::string pFstr = pathstr + ".implicit.cst.smb.Fp.txt";
  std::string Sstr = pathstr + ".implicit.cst.smb.S.txt";
  std::string KSstr = pathstr + ".implicit.cst.smb.KS.txt";
  std::string sigmastr = pathstr + ".implicit.cst.smb.sigma.txt";
  std::string pcstr = pathstr + ".implicit.cst.smb.pc.txt";
  std::string matstr = pathstr + ".implicit.cst.smb.mat.txt";
  
  // Time integrator
  double Dt = 0.01, InitialTime=0, FinalTime = 13200 * Dt, CurrentTime=0;
  size_t TimeStep = 0;
  TimeIntegrationDataManager TimeIntegrationData( InitialTime, FinalTime, Dt, CurrentTime, TimeStep );
  
  // Parameters and Model definition
  bool usingSmoothMacauleyBrackets = true;
  KMS_IJSS2017_Parameters matmodel( 1.0, 1.1, 0.15, 0.0005, 0.62, 0.37, 77.22, 13.01, 15, 0.01, 0.2, 5.8, 30, 60, 0.063, 0.008, 2, 1, 290, usingSmoothMacauleyBrackets ) ;
  KMS_IJSS2017_IO IO(logstr, Fstr, pFstr, Sstr, KSstr, sigmastr, pcstr, matstr, PrintEveryNSteps);
  MechanicalModels_Test_Functions<spacedim> KMS_IJSS2017TestFunctions;
  
  KMS_IJSS2017_Implicit_BE_Staggered<spacedim> ImplicitModelInstance( &matmodel, &IO, &TimeIntegrationData );
  
  
  // Consistent Tangent test
  bool Verbose = false;
  
  t = clock();
  IO.log()  << "\n  ImplicitModelInstance.IntegratorTest for CompactionAndShear1 running. \n";
  ImplicitModelInstance.ConsistentTangentTest( KMS_IJSS2017TestFunctions.IsotropicCompaction, Verbose );
  t = clock() - t;
  IO.log()  << "\n  ImplicitModelInstance.IntegratorTest for CompactionAndShear1 using Smooth Macauley Brackets terminated in " << ((float)t)/CLOCKS_PER_SEC << "seconds. \n";
  
}

int main(int argc, const char * argv[])
{  
  system("mkdir -p out/KMS_IJSS2017");
  
  std::ofstream logfile( "./ttl-learning.log", std::ios::out );
  logfile << "\n\n" << "ttl-learning started properly. \n\n" << std::flush;
    
//  logfile << "\n KMS_IJSS2017 Explicit test running " << std::flush;
//  KMS_IJSS2017_Explicit_test();
//  logfile << "/ completed \n" << std::flush;
  
  logfile << "\n KMS_IJSS2017 Implicit test running " << std::flush;
  KMS_IJSS2017_Implicit_test();
  logfile << "/ completed \n\n" << std::flush;
  
//  logfile << "\n KMS_IJSS2017 Implicit test running " << std::flush;
//  KMS_IJSS2017_Implicit_ConsistentTangent_test();
//  logfile << "/ completed \n\n" << std::flush;

  logfile << "ttl-learning terminated correctly. \n" << std::flush;
  
  logfile.close();
  
  return 0;
}










