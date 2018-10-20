#ifndef H__H__HYPERELASTICITY__H__H
#define H__H__HYPERELASTICITY__H__H

#include "material_properties.h"


typedef void (*deviatoric_part)(double *C_in, MATERIAL_ELASTICITY const *mat, double *S);
typedef void (*volume_part)(double *dudj, double J);

class StrainEnergyDensityFunction
{
  public:
    deviatoric_part compute_potential_dev;  
    deviatoric_part compute_PK2_dev;
    deviatoric_part compute_tangent_dev;
    deviatoric_part compute_d3W_dC3_dev;  
    volume_part compute_u;  
    volume_part compute_dudj;
    volume_part compute_d2udj2;
    volume_part compute_d3udj3;    
    
    /// constructor
    StrainEnergyDensityFunction(){ set_null(); };
    
    /// set strain energy density functions
    int set_material(const MATERIAL_ELASTICITY *mat);

    /// compute PKII stress and elasticity tensor (4th order) and update references (S_out, and L_out)
    /// rather than member S and L 
    /// if update_stiffness is false, there will be no computing elasticity tensor
    /// s.t NULL is valid for L_out
    int update_elasticity(const MATERIAL_ELASTICITY *mat,
                          double *S_out,
                          double *L_out,
                          double *Fe, 
                          const bool update_stiffness);
                          
    /// compute 3rd derivative of strain energy function w.r.t eC        
    int compute_d3W_dC3(const MATERIAL_ELASTICITY *mat,
                        double *eF, 
                        double *K_out);
    /// set null
    void set_null(){
      compute_potential_dev = NULL;
      compute_PK2_dev       = NULL;
      compute_tangent_dev   = NULL;
      compute_d3W_dC3_dev   = NULL;
      compute_u             = NULL;
      compute_dudj          = NULL;
      compute_d2udj2        = NULL;
      compute_d3udj3        = NULL;
    };
    
    ~StrainEnergyDensityFunction(){ set_null(); };
};

/// comute Von Mises stress
double compute_Von_Mises(double *T_in);

/// compute Cauchy stress from PKII in elastic configurateion
void compute_Cauchy_stress(double *PKII,
                           double *sigma_out, 
                           double *eF);

template<class Mat, class Sedf> class GcmElasticity
{
  public:
    // if ture, S and L are created
    bool is_SL_created;
    
    // memory for PKII and 4th order elasticity tensor
    double *L, *S;
    
    // address for material object
    Mat *mat;
    
    // define strain energy density function
    Sedf sedf;
    
    // if ture compute 4th order elasticity tensor
    bool compute_stiffness;
    
    GcmElasticity(){
      L = S = NULL;
      mat = NULL;
      is_SL_created     = false;
      compute_stiffness = false;
    };
    
    GcmElasticity(Mat *mat, 
                  const bool cp_stiff = false){
      construct_elasticity(mat, cp_stiff);
    };       
    
    ~GcmElasticity(){
      if(is_SL_created){
        delete [] S;
        if(compute_stiffness)
          delete [] L;
        is_SL_created = false;
      }
    }
    
    /// allocate memories and assign material object
    /// 
    /// \param[in] mat      material object
    /// \param[in] cp_stiff if true compute stiffness (4th order elasticity tensor)
    /// \param[in] cp_stiff if false no compute stiffness
    /// \return non-zero on internal error
    int construct_elasticity(Mat *mat, const bool cp_stiff){
      int err = 0;
      this->is_SL_created = true;
    
      this->S = new double [9];
      if(cp_stiff)
        this->L = new double [81];
  
      this->mat  = mat;
      this->compute_stiffness = cp_stiff;
      this->set_SEDF();
      return err;
    };
      
    int set_SEDF(void){
      return sedf.set_material(this->mat);
    }

    //int construct_elasticity(Mat *mat, const bool cp_stiff = false);

    /// compute PKII stress and elasticity tensor (4th order) and update references 
    /// (S_out, and L_out) rather than member S and L
    int update_elasticity(double *S_out,
                          double *L_out,
                          double *Fe, 
                          const bool update_stiffness){
      return sedf.update_elasticity(mat, S_out, L_out, Fe, update_stiffness);
    }

    /// compute PKII stress and elasticity tensor (4th order)
    /// and update member S and L 
    /// if update_stiffness is false, there will be no computing elasticity tensor
    /// 
    /// \param[in]  Fe    elastic deformation gradient
    /// \param[in]  update_stiffness if true, compute 4th order elasticity Tensor
    ///                                 false no compute elasticity Tensor
    /// \return non-zero on internal error
    int update_elasticity(double *Fe, 
                          const bool update_stiffness){
      return this->update_elasticity(this->S, this->L, Fe, update_stiffness);
    };

    /// compute PKII stress and 
    /// if member compute_stiffness is ture compute elasticity tensor (4th order)
    /// 
    /// \param[in]  Fe elastic deformation gradient
    /// \return non-zero on internal error
    int update_elasticity(double *Fe){
      return this->update_elasticity(this->S, this->L, Fe, this->compute_stiffness);
    };
    
    /// compute deviatoric part of potential energy
    ///
    /// \param[in]  C_in right Caucy Green tensor (2nd order)
    /// \param[out] P    computed potential energy
    void compute_potential_dev(double *C_in, double *P){
      this->sedf.compute_potential_dev(C_in, this->mat, P);
    };

    /// compute deviatoric part of PKII
    ///
    /// \param[in]  C_in  right Caucy Green tensor (2nd order)
    /// \param[out] S_out computed PKII (2nd order)   
    void compute_PK2_dev(double *C_in, double *S_out){
      this->sedf.compute_PK2_dev(C_in, this->mat, S_out);
    };

    /// compute deviatoric part of elasticity tensor
    ///
    /// \param[in]  C_in  right Caucy Green tensor (2nd order)
    /// \param[out] L_out computed elasticity tensor (4th order)
    void compute_tangent_dev(double *C_in, double *L_out){
      this->sedf.compute_tangent_dev(C_in, this->mat, L_out);
    };
    
    /// compute 3rd derivative of strain energy function (deviatoric part) w.r.t eC
    ///
    /// \param[in]  C_in right Caucy Green tensor (2nd order)
    /// \param[out] I    computed d3WdC3 tensor (6th order)    
    void compute_d3W_dC3_dev(double *C_in, double *K){
      this->sedf.compute_d3W_dC3_dev(C_in, this->mat, K);
    };

    /// compute volumetric part of potential energy
    ///
    /// \param[out] u computed potential energy    
    /// \param[in]  J determinant of elastic deformation gradient tensor
    void compute_u(double *u, const double J){
      this->sedf.compute_u(u, J);
    };

    /// compute 1st derivative of volumetric part of potential energy w.r.t J
    ///
    /// \param[out] du computed derivative of volumetric part of potential energy
    /// \param[in]  J  determinant of elastic deformation gradient tensor    
    void compute_dudj(double *du, const double J){
      this->sedf.compute_dudj(du, J);
    };

    /// compute 2nd derivative of volumetric part of potential energy w.r.t J
    ///
    /// \param[out] ddu computed derivative of volumetric part of potential energy
    /// \param[in]  J   determinant of elastic deformation gradient tensor 
    void compute_d2udj2(double *ddu, const double J){
      this->sedf.compute_d2udj2(ddu, J);
    };

    /// compute 3rd derivative of volumetric part of potential energy w.r.t J
    ///
    /// \param[out] d3u computed derivative of volumetric part of potential energy
    /// \param[in]  J   determinant of elastic deformation gradient tensor 
    void compute_d3udj3(double *d3u, const double J){
      this->sedf.compute_d3udj3(d3u, J);
    };

    /// comute Von Mises of PKII
    /// internal S(PKII) will be used
    ///
    /// return computed PKII^eff
    double compute_PK2_eff(void){
      return compute_Von_Mises(this->S);
    };   

    /// comute Von Mises of Cauchy stress
    /// internal S(PKII) will be used to compute Cauchy stress
    ///
    /// \param[out] sigma_eff computed sigma^eff
    /// \param[in]  eF        elastic deformation gradient
    /// \return non-zero on internal error
    void compute_Cauchy_eff(double *sigma_eff, double *eF){
      double sigma[9];  
      this->compute_Cauchy(sigma, eF);
      *sigma_eff = compute_Von_Mises(sigma);
    };
      
    /// comute Cauchy stress
    /// internal S(PKII) will be used to compute Cauchy stress
    ///
    /// \param[out] sigma computed Chauchy stress
    /// \param[in]  eF    elastic deformation gradient
    void compute_Cauchy(double *sigma_out, double *eF){
      compute_Cauchy_stress(this->S, sigma_out, eF);
    };
    
    int compute_d3W_dC3(double *eF, double *K_out){
      return sedf.compute_d3W_dC3(mat, eF, K_out);
    };
};

class HyperElasticity: public GcmElasticity<MATERIAL_ELASTICITY,StrainEnergyDensityFunction>
{};


#endif 
