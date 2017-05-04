//
//  KMS-IJSS2017-parameters.h
//  ttl-learning
//
//  Created by alberto salvadori on 5/3/17.
//  Copyright © 2017 alberto salvadori. All rights reserved.
//

#ifndef KMS_IJSS2017_parameters_h
#define KMS_IJSS2017_parameters_h


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



#endif /* KMS_IJSS2017_parameters_h */
