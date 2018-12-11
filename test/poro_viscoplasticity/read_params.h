#ifndef __H__GCM_PVP_TEST__READ_PARAMETERS_H__
#define __H__GCM_PVP_TEST__READ_PARAMETERS_H__

enum param_names {
  PARAM_yf_M,       // Yield function parameters
  PARAM_yf_alpha,   //   :
  PARAM_flr_m_d,    // Flow rule parameters
  PARAM_flr_m_v,    // Flow rule parameters  
  PARAM_flr_gamma0, //   :
  PARAM_hr_a1,      // Hardening rule parameters
  PARAM_hr_a2,      //   :
  PARAM_hr_Lambda1, //   :
  PARAM_hr_Lambda2, //   :
  PARAM_c_inf,      // Cohesion rule parameters
  PARAM_c_Gamma,    //   :
  PARAM_d_B,        // Transition rule parameters
  PARAM_d_pcb,      //   :
  PARAM_d_m,        //   :
  PARAM_mu_0,       // Shear modulus parameters
  PARAM_mu_1,       //   :
  PARAM_K_p0,       // Bulk modulus parameters
  PARAM_K_kappa,    //   :
  PARAM_pl_n,       // Power law exponent
  PARAM_cf_g0,      // Compaction function parameters
  PARAM_cf_pcinf,   //   :
  PARAM_NO
};


/// read poro-visco-plasticity parameters from a file input formatted as
/// M alpha m_d m_v gamma0 a1 a2 Lambda1 Lambda2 inf Gamma B pcb d_m mu_0 mu_1 p0 kappa n g0 pc_inf
/// 
/// \param[in]      *fp      file pointer to be read
/// \param[in, out] &mat_pvp material property object for porovisco-plasticity
void read_pvp_material_properties(FILE *fp,
                                  MaterialPoroViscoPlasticity &mat_pvp){
  // read material parameters 
  double param[PARAM_NO];
  for(int ia=0; ia<PARAM_NO; ia++)
    fscanf(fp, "%lf", param+ia);
    
  set_properties_poro_visco_plasticity(&mat_pvp,
                                       param[PARAM_yf_M],
                                       param[PARAM_yf_alpha],
                                       param[PARAM_flr_m_d],
                                       param[PARAM_flr_m_v],
                                       param[PARAM_flr_gamma0],
                                       param[PARAM_hr_a1],
                                       param[PARAM_hr_a2],
                                       param[PARAM_hr_Lambda1],
                                       param[PARAM_hr_Lambda2],
                                       param[PARAM_c_inf],
                                       param[PARAM_c_Gamma],
                                       param[PARAM_d_B],
                                       param[PARAM_d_pcb],
                                       param[PARAM_d_m],
                                       param[PARAM_mu_0],
                                       param[PARAM_mu_1],
                                       param[PARAM_K_p0],
                                       param[PARAM_K_kappa],
                                       param[PARAM_pl_n],
                                       param[PARAM_cf_g0],
                                       param[PARAM_cf_pcinf]);  
}

#endif