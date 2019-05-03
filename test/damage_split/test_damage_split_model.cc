#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "continuum_damage_model.h"
#include <sys/time.h>

#include <stdio.h>
#include <math.h>
#define PI 3.141592653589793

#define F_ROW_MAX 300
#define F_COLUMN_MAX 10
#define DIM_3   3
#define DIM_3x3 9

typedef int (*DEFORMATION_GRADIENT) (Tensor<2> &F, double t); 

typedef struct {
  double E;
  double nu;
  double p1;
  double p2;
  double Yin;
  double mu;
  double w_max;
  double alpha_dev;
  double beta_dev_p;
  double beta_dev_m;  
  double alpha_vol;
  double beta_vol_p;
  double beta_vol_m;  
} MAT_PROP;

typedef struct {
  double dt;
  int stepno;
  char file_out[1024];
  DEFORMATION_GRADIENT F_of_t;
} SIM_PARAMS;

int F_of_t_tension_const_J(Tensor<2> &F, double t)
{
  int err = 0;
  double d = 0.0002;
  double F11 = 1.0 + d*t;
  F.data[0] = F11;
  F.data[4] = F.data[8] = 1.0/sqrt(F11);
    
  F.data[1] = F.data[2] = F.data[3] = 0.0;
  F.data[5] = F.data[6] = F.data[7] = 0.0;
  
  return err;
};

int F_of_t_extension(Tensor<2> &F, double t)
{
  int err = 0;
  double d = 0.0002;
  double F11 = 1.0 + d*t;
  F.data[0] = F.data[4] = F.data[8] = F11;    
  F.data[1] = F.data[2] = F.data[3] = 0.0;
  F.data[5] = F.data[6] = F.data[7] = 0.0;
  
  return err;
};

int F_of_t_mixed(Tensor<2> &F, double t)
{
  int err = 0;
  double d = 0.0005;
  double t1 = 25.0;
  double t2 = 50.0;
  double t3 = 75.0;

  F.data[0] = F.data[4] = F.data[8] = 1.0;    
  F.data[1] = F.data[2] = F.data[3] = 0.0;
  F.data[5] = F.data[6] = F.data[7] = 0.0;
      
  if(t<=t1)
  {
    double F11 = d*t;
    F.data[0] = 1.0 + F11;
  }
  if(t1<t && t<=t2)
  {
    double F11 = d*t1;
    double F21 = d*(t-t1);
    F.data[0] = 1.0 + F11;
    F.data[3] = F21; 
  }
  if(t2<t && t<=t3)
  {
    double F11 = d*t1;
    double F21 = d*(t2-t1) - d*(t-t2);
    F.data[0] = 1.0 + F11;
    F.data[3] = F21; 
  }
    if(t3<t)
  {
    double F11 = d*t1;
    double F22 = d*(t-t3);
    F.data[0] = 1.0 + F11;
    F.data[4] = 1.0 + F22;;    
  }  
  
  return err;
};

int F_of_t_tension_compression(Tensor<2> &F, double t)
{
  int err = 0;
  double m = 0.02;
  double tmax = 100;
  double T = 2;
  double F11 = m*sin(t/tmax*2.0*T*PI);
  F.data[0] = 1.0 + F11;
  F.data[4] = F.data[8] = 1.0 + 0.25*F11;
    
  F.data[1] = F.data[2] = F.data[3] = 0.0;
  F.data[5] = F.data[6] = F.data[7] = 0.0;
  
  return err;
};

int read_deformation_gradient(Tensor<2> &F)
{
  int err = 0;
  FILE *fp = fopen("deformation_gradient.txt", "r");

  for(int a = 0; a<F_ROW_MAX; a++)
  {
    for(int b=0; b<F_COLUMN_MAX; b++)
      fscanf(fp, "%lf", (F.data) + a*F_COLUMN_MAX + b);
  }
  fclose(fp);
  return err;
};
                               
int print_output_header(FILE *out)
{
  int vno = 0;
  fprintf(out,"%% %s[%d] %s[%d] ", "t", vno, "Jnp1", vno+1);
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ", 
               "F11", vno+2, "F12", vno+3, "F13",vno+4,
               "F21", vno+5, "F22", vno+6, "F23",vno+7,
               "F31", vno+8, "F32", vno+9, "F33",vno+10);
  vno += 10;               
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] ",
               "dw", vno+1,"vw", vno+2,"dX", vno+3,"vX", vno+4);
  vno += 4;               
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ",
               "sigma0_11", vno+1,"sigma0_22", vno+2,"sigma0_33", vno+3,
               "sigma0_23", vno+4,"sigma0_31", vno+5,"sigma0_12", vno+6,"sigma_eff_0",vno+7);
  vno += 7;               
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ",
               "S0_11", vno+1, "S022", vno+2, "S033",vno+3,
               "S0_23", vno+4, "S031", vno+5, "S012",vno+6);
  vno += 6;
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ",
               "sigma11", vno+1, "sigma22", vno+2, "sigma33",vno+3,
               "sigma23", vno+4, "sigma31", vno+5, "sigma12",vno+6,"sigma_eff",vno+7);
  vno += 7;
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ",
               "S11", vno+1, "S22", vno+2, "S33",vno+3,
               "S23", vno+4, "S31", vno+5, "S12",vno+6);
  vno += 6;
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ",
               "E11", vno+1, "E22", vno+2, "E33",vno+3,
               "E23", vno+4, "E31", vno+5, "E12",vno+6,"Eeff",vno+7);
  vno += 7;
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d]\n",
               "e11", vno+1, "e22", vno+2, "e33",vno+3,
               "e23", vno+4, "e31", vno+5, "e12",vno+6,"eeff",vno+7);  
  return 0;
};

int print_outputs(FILE *out, double t,
                             double J,
                             double dw, double vw, double dX, double vX,
                             double sigma_eff_0, double sigma_eff, double Eeff, double eeff,
                             double *F,  double *sigma0,
                             double *S0, double *sigma,
                             double *S,  double *E,
                             double *e)
{
  fprintf(out,"%e %e ", t, J);
  fprintf(out,"%e %e %e %e %e %e %e %e %e ", 
               F[0],F[1],F[2],F[3],F[4],F[5],F[6],F[7],F[8]);
  fprintf(out,"%e %e %e %e ", dw,vw,dX,vX);
  fprintf(out,"%e %e %e %e %e %e %e ",
               sigma0[0],sigma0[4],sigma0[8],sigma0[5],sigma0[6],sigma0[1], sigma_eff_0);
  fprintf(out,"%e %e %e %e %e %e ",
               S0[0],S0[4],S0[8],S0[5],S0[6],S0[1]);
  fprintf(out,"%e %e %e %e %e %e %e ",
               sigma[0],sigma[4],sigma[8],sigma[5],sigma[6],sigma[1], sigma_eff);
  fprintf(out,"%e %e %e %e %e %e ",
               S[0],S[4],S[8],S[5],S[6],S[1]);
  fprintf(out,"%e %e %e %e %e %e %e ",
               E[0],E[4],E[8],E[5],E[6],E[1],Eeff);
  fprintf(out,"%e %e %e %e %e %e %e\n",
               e[0],e[4],e[8],e[5],e[6],e[1],eeff);  
  return 0;
};

int test_damage_model(MAT_PROP *mat, SIM_PARAMS *sim)
{
  int err = 0;

  // ----------------------------------------------------------
  // construct material property structures and initialize them
  // ----------------------------------------------------------
  MATERIAL_ELASTICITY        mat_e; // for elasticity
  MATERIAL_CONTINUUM_DAMAGE  mat_d; // for damage
  err += set_properties_using_E_and_nu(&mat_e,mat->E,mat->nu);

  err += set_damage_parameters(&mat_d, mat->p1, 
                                       mat->p2, 
                                       mat->Yin, 
                                       mat->mu, 
                                       mat->w_max);
  
  // --------------------------------------------------
  // construct and initialize elasticity handle
  // : provides computing manner of elasticity response
  // --------------------------------------------------
  bool construct_ET = false; // computing elasticity tensor
  HyperElasticity elast;
  elast.construct_elasticity(&mat_e, construct_ET);
    
  // -----------------------
  // construct GL(3) tensors
  // -----------------------
  Tensor<2> Fnp1,E,sigma0,S0,sigma,e,FI,eye;
  eye = ttl::identity(i,j);
    
  Fnp1 = eye(i,j);

  TensorA<2> S(elast.S);

  // --------------
  // do simulations
  // --------------
  
  double Jnp1;
  double dt = sim->dt;  
  FILE *out = fopen(sim->file_out, "w");
  err += print_output_header(out);
 
  double w, X, H, wn, Xn;
  w = wn = X = Xn = H = 0.0;
  int is_it_damaged = 0;

  for(int a = 0; a<=sim->stepno; a++)
  {
    double t = dt*a;
    
    err += (sim->F_of_t)(Fnp1, t);
    
    // compute strains
    E = 0.5*(Fnp1(k,i)*Fnp1(k,j) - eye(i,j));

    double Eeff = E(i,j)*E(i,j);
    Eeff = sqrt(Eeff);
    
    err += inv(Fnp1,FI);
    e = FI(k,i)*E(k,l)*FI(j,l);
    double eeff = e(i,j)*e(i,j);
    eeff = sqrt(eeff);

        
    //compute undamaged values
    double sigma_eff_0 = 0.0;
    
    err += elast.update_elasticity(Fnp1.data, construct_ET);
    elast.compute_Cauchy_eff(&sigma_eff_0, Fnp1.data);
    elast.compute_Cauchy(sigma0.data, Fnp1.data);
    
    S0 = S(i,j);
        
    //compute damaged values
    double sigma_eff = 0.0;
    Jnp1 = ttl::det(Fnp1);

    err += continuum_damage_integration_alg(&mat_d,&elast,
                                            &w,&X,&H,&is_it_damaged,
                                            wn,Xn,dt,Fnp1.data);
                                                                                 
    err += update_damage_elasticity(&mat_d,&elast,w,is_it_damaged,H,
                                    dt,Fnp1.data, construct_ET);
                                            
    elast.compute_Cauchy_eff(&sigma_eff,Fnp1.data);
    elast.compute_Cauchy(sigma.data,Fnp1.data);
    
    err += print_outputs(out,t,Jnp1,w,0.0,X,0.0,
                                  sigma_eff_0,sigma_eff,Eeff,eeff,
                                  Fnp1.data, sigma0.data,
                                    S0.data,  sigma.data,
                                     S.data,      E.data,
                                     e.data);
    wn = w;
    Xn = X;
    is_it_damaged = 0;
  }
    
  return err;    
}

int test_split_damage_model(MAT_PROP *mat, SIM_PARAMS *sim)
{
  int err = 0;

  // ----------------------------------------------------------
  // construct material property structures and initialize them
  // ----------------------------------------------------------
  MATERIAL_ELASTICITY        mat_e; // for elasticity
  MATERIAL_CONTINUUM_DAMAGE  mat_d; // for damage
  err += set_properties_using_E_and_nu(&mat_e,mat->E,mat->nu);


  err += set_split_damage_parameters(&mat_d, mat->p1, 
                                             mat->p2, 
                                             mat->Yin, 
                                             mat->mu, 
                                             mat->w_max,
                                             mat->alpha_dev,
                                             mat->beta_dev_p,
                                             mat->beta_dev_m,                                             
                                             mat->alpha_vol,
                                             mat->beta_vol_p,
                                             mat->beta_vol_m);                                             

  //err += print_material_property_elasticity(&mat_e); // just printing
  
  // --------------------------------------------------
  // construct and initialize elasticity handle
  // : provides computing manner of elasticity response
  // --------------------------------------------------
  bool construct_ET = true; // computing elasticity tensor
  HyperElasticity elast;
  elast.construct_elasticity(&mat_e, construct_ET);
    
  // -----------------------
  // construct GL(3) tensors
  // -----------------------
  Tensor<2> F,Fnp1,E,sigma0,S0,sigma,e,FI,eye;
  
  eye = ttl::identity(i,j);
  TensorA<2> S(elast.S);;

  // --------------
  // do simulations
  // --------------
  
  double Jnp1;
  double dt = sim->dt;  

  double dw,vw,dX,vX,dH,vH;
  int is_it_damaged_d, is_it_damaged_v;
  double dwn,vwn,dXn,vXn; 
  
  dw=vw=dX=vX=dH=vH=dwn=vwn=dXn=vXn=0.0;
  is_it_damaged_d = is_it_damaged_v = 0;


  FILE *out = fopen(sim->file_out, "w");
  err += print_output_header(out);

  for(int a = 0; a<=sim->stepno; a++)
  {
    double t = dt*a;
    
    err += (sim->F_of_t)(Fnp1, t);
    
    // compute strains
    // compute strains
    E = 0.5*(Fnp1(k,i)*Fnp1(k,j) - eye(i,j));

    double Eeff = E(i,j)*E(i,j);
    Eeff = sqrt(Eeff);
    
    err += inv(Fnp1,FI);
    e = FI(k,i)*E(k,l)*FI(j,l);
    double eeff = e(i,j)*e(i,j);
    eeff = sqrt(eeff);
        
    //compute undamaged values
    double sigma_eff_0 = 0.0;
    
    err += elast.update_elasticity(Fnp1.data, construct_ET);
    elast.compute_Cauchy_eff(&sigma_eff_0, Fnp1.data);
    elast.compute_Cauchy(sigma0.data, Fnp1.data);
    S0 = S(i,j);
        
    //compute damaged values
    double sigma_eff = 0.0;
    Jnp1 = ttl::det(Fnp1);
    
    err += continuum_split_damage_integration_alg(&mat_d,&elast,
                                     &dw,&vw,&dX,&vX,&dH,&vH,
                                     &is_it_damaged_d,&is_it_damaged_v,                                     
                                     dwn,vwn,dXn,vXn,dt,Fnp1.data);    
    
    err += update_split_damage_elasticity(&mat_d,&elast,dw,vw,dH,vH,
                                          is_it_damaged_d,is_it_damaged_v,
                                          dt, Fnp1.data, construct_ET);
                                          
    elast.compute_Cauchy_eff(&sigma_eff,Fnp1.data);
    elast.compute_Cauchy(sigma.data,Fnp1.data);
    
    err += print_outputs(out,t,Jnp1,dw,vw,dX,vX,
                                  sigma_eff_0,sigma_eff,Eeff,eeff,
                                  Fnp1.data, sigma0.data,
                                    S0.data,  sigma.data,
                                     S.data,      E.data,
                                     e.data);

    dwn = dw;
    vwn = vw;
    dXn = dX;
    vXn = vX;
    
    is_it_damaged_d = 0;
    is_it_damaged_v = 0;
  } 
    
  return err;
}
int run_simulations(void)
{
  int err = 0;
  
  // set material properties
  MAT_PROP mat;  
  
  mat.E       = 9.95e+3;
  mat.nu      = 0.34;
  mat.p1      = 8.0;
  mat.p2      = 2.5;
  mat.Yin     = 0.15;
  mat.mu      = 0.1;
  mat.w_max   = 0.98;
  mat.alpha_dev  = 1.0;
  mat.beta_dev_p = 0.0;
  mat.beta_dev_m = 0.0;  
  mat.alpha_vol  = 0.0;
  mat.beta_vol_p = 1.0;
  mat.beta_vol_m = 1.0;  

  // set simulation parameters
  SIM_PARAMS sim;
  sim.dt     = 1.0;
  sim.stepno = 100;  
  sim.F_of_t = F_of_t_tension_const_J;  
  sprintf(sim.file_out, "tension_const_J_mu_0p1.txt");
  err += test_split_damage_model(&mat, &sim);
  
  mat.mu = 100.0;
  sprintf(sim.file_out, "tension_const_J_mu_100p0.txt");
  err += test_split_damage_model(&mat, &sim);
  
  mat.mu = 0.01;
  sprintf(sim.file_out, "tension_const_J_mu_0p01.txt");
  err += test_split_damage_model(&mat, &sim);

  mat.mu = 0.1;
  
  mat.alpha_vol = 0.3;
  sprintf(sim.file_out, "tension_const_J_alpha_vol_0p3.txt");
  err += test_split_damage_model(&mat, &sim);  

  mat.alpha_vol = 0.7;
  sprintf(sim.file_out, "tension_const_J_alpha_vol_0p7.txt");
  err += test_split_damage_model(&mat, &sim);
  
  mat.mu = 0.1;
  mat.alpha_vol = 0.0;
  sim.F_of_t = F_of_t_extension;

  mat.beta_dev_p = 0.0;
  mat.beta_dev_m = 0.0;  
  sprintf(sim.file_out, "expansion_beta_dev_0p0.txt");
  err += test_split_damage_model(&mat, &sim);  

  mat.beta_dev_p = 0.3;
  mat.beta_dev_m = 0.3;  
  sprintf(sim.file_out, "expansion_beta_dev_0p3.txt");
  err += test_split_damage_model(&mat, &sim);  

  mat.beta_dev_p = 0.7;
  mat.beta_dev_m = 0.7;  
  sprintf(sim.file_out, "expansion_beta_dev_0p7.txt");
  err += test_split_damage_model(&mat, &sim);

  mat.mu = 0.1;
  mat.alpha_vol = 0.0;
  mat.beta_dev_p = 0.0; 
  mat.beta_dev_m = 0.0;   
  sim.F_of_t = F_of_t_mixed;

  sprintf(sim.file_out, "tension_mu_0p1.txt");
  err += test_split_damage_model(&mat, &sim);  

  mat.mu = 100.0;
  sprintf(sim.file_out, "tension_mu_100p0.txt");
  err += test_split_damage_model(&mat, &sim);

  mat.mu = 0.01;
  sprintf(sim.file_out, "tension_mu_0p01.txt");
  err += test_split_damage_model(&mat, &sim);
  
  mat.mu = 0.1;
  mat.alpha_vol = 0.2;
  mat.beta_dev_p = 0.0;
  mat.beta_dev_m = 0.0;
  sprintf(sim.file_out, "tension_alpha_vol_0p2.txt");
  err += test_split_damage_model(&mat, &sim);  

  mat.alpha_vol = 0.0;
  mat.beta_dev_p = 0.2;
  mat.beta_dev_m = 0.2;
  sprintf(sim.file_out, "tension_beta_dev_0p2.txt");
  err += test_split_damage_model(&mat, &sim);
  
  mat.alpha_vol = 0.2;
  mat.beta_dev_p = 0.2;
  mat.beta_dev_m = 0.2;
  sprintf(sim.file_out, "tension_alpha_vol_0p2_beta_dev_0p2.txt");
  err += test_split_damage_model(&mat, &sim);
 
  mat.alpha_vol = 1.0;
  mat.beta_dev_p = 1.0;
  mat.beta_dev_m = 1.0;  
  sprintf(sim.file_out, "tension_alpha_vol_1p0_beta_dev_1p0.txt");
  err += test_split_damage_model(&mat, &sim);

  sprintf(sim.file_out, "tension_original_damage.txt");
  err += test_damage_model(&mat, &sim);
  
  sim.F_of_t = F_of_t_tension_compression;
  mat.alpha_vol = 0.2;
  mat.beta_dev_p = 0.2;
  mat.beta_dev_m = 0.2;
  sprintf(sim.file_out, "tension_compression_beta_dev_0p2_beta_vol_1p0.txt");
  err += test_split_damage_model(&mat, &sim);
  
  mat.beta_vol_p = 1.0;
  mat.beta_vol_m = 1.0;  
  mat.beta_dev_p = 0.0;
  mat.beta_dev_m = 0.0;  
  sprintf(sim.file_out, "tension_compression_beta_dev_0p0_beta_vol_m1p0.txt");
  err += test_split_damage_model(&mat, &sim);
  
  mat.beta_vol_p = 1.0;
  mat.beta_vol_m = 1.0;  
  mat.beta_dev_p = 0.2;
  mat.beta_dev_m = 0.2;

  sprintf(sim.file_out, "tension_compression_beta_dev_m0p2_beta_vol_m1p0.txt");
  err += test_split_damage_model(&mat, &sim);      
      
  return err;
}


int main(int argc,char *argv[])
{
  int err = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  err += run_simulations();    
  
  gettimeofday(&end, NULL);
  double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0 
  + (double)(end.tv_sec - start.tv_sec);
  printf ("Total time: %.4lf s\n", diff);
  
  return err;
}

