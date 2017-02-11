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

typedef int (*DEFORMATION_GRADIENT) (Matrix(double) *F, double t); 

typedef struct {
  double E;
  double nu;
  double p1;
  double p2;
  double Yin;
  double mu;
  double w_max;
  double alpha_dev;
  double beta_dev;
  double alpha_vol;
  double beta_vol;  
} MAT_PROP;

typedef struct {
  double dt;
  int stepno;
  char file_out[1024];
  DEFORMATION_GRADIENT F_of_t;
} SIM_PARAMS;

int F_of_t_tension_const_J(Matrix(double) *F, double t)
{
  int err = 0;
  double d = 0.0002;
  double F11 = 1.0 + d*t;
  F->m_pdata[0] = F11;
  F->m_pdata[4] = F->m_pdata[8] = 1.0/sqrt(F11);
    
  F->m_pdata[1] = F->m_pdata[2] = F->m_pdata[3] = 0.0;
  F->m_pdata[5] = F->m_pdata[6] = F->m_pdata[7] = 0.0;
  
  return err;
};

int F_of_t_extension(Matrix(double) *F, double t)
{
  int err = 0;
  double d = 0.0002;
  double F11 = 1.0 + d*t;
  F->m_pdata[0] = F->m_pdata[4] = F->m_pdata[8] = F11;    
  F->m_pdata[1] = F->m_pdata[2] = F->m_pdata[3] = 0.0;
  F->m_pdata[5] = F->m_pdata[6] = F->m_pdata[7] = 0.0;
  
  return err;
};

int F_of_t_mixed(Matrix(double) *F, double t)
{
  int err = 0;
  double d = 0.0005;
  double t1 = 25.0;
  double t2 = 50.0;
  double t3 = 75.0;

  F->m_pdata[0] = F->m_pdata[4] = F->m_pdata[8] = 1.0;    
  F->m_pdata[1] = F->m_pdata[2] = F->m_pdata[3] = 0.0;
  F->m_pdata[5] = F->m_pdata[6] = F->m_pdata[7] = 0.0;
      
  if(t<=t1)
  {
    double F11 = d*t;
    F->m_pdata[0] = 1.0 + F11;
  }
  if(t1<t && t<=t2)
  {
    double F11 = d*t1;
    double F21 = d*(t-t1);
    F->m_pdata[0] = 1.0 + F11;
    F->m_pdata[3] = F21; 
  }
  if(t2<t && t<=t3)
  {
    double F11 = d*t1;
    double F21 = d*(t2-t1) - d*(t-t2);
    F->m_pdata[0] = 1.0 + F11;
    F->m_pdata[3] = F21; 
  }
    if(t3<t)
  {
    double F11 = d*t1;
    double F22 = d*(t-t3);
    F->m_pdata[0] = 1.0 + F11;
    F->m_pdata[4] = 1.0 + F22;;    
  }  
  
  return err;
};

int F_of_t_tension_compression(Matrix(double) *F, double t)
{
  int err = 0;
  double m = 0.02;
  double tmax = 100;
  double T = 2;
  double F11 = m*sin(t/tmax*2.0*T*PI);
  F->m_pdata[0] = 1.0 + F11;
  F->m_pdata[4] = F->m_pdata[8] = 1.0 + 0.25*F11;
    
  F->m_pdata[1] = F->m_pdata[2] = F->m_pdata[3] = 0.0;
  F->m_pdata[5] = F->m_pdata[6] = F->m_pdata[7] = 0.0;
  
  return err;
};
  
int read_deformation_gradient(Matrix(double) *F)
{
  int err = 0;
  FILE *fp = fopen("deformation_gradient.txt", "r");

  for(int a = 0; a<F_ROW_MAX; a++)
  {
    for(int b=0; b<F_COLUMN_MAX; b++)
      fscanf(fp, "%lf", (F->m_pdata) + a*F_COLUMN_MAX + b);
  }
  fclose(fp);
  return err;
};
                               
int print_output_header(FILE *out)
{
  int vno = 0;
  fprintf(out,"%% %s[%d] %s[%d] ", "t", vno++, "Jnp1", vno++);
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ", 
               "F11", vno++, "F12", vno++, "F13",vno++,
               "F21", vno++, "F22", vno++, "F23",vno++,
               "F31", vno++, "F32", vno++, "F33",vno++);
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] ",
               "dw", vno++,"vw", vno++,"dX", vno++,"vX", vno++);
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ",
               "sigma0_11", vno++,"sigma0_22", vno++,"sigma0_33", vno++,
               "sigma0_23", vno++,"sigma0_31", vno++,"sigma0_12", vno++,"sigma_eff_0",vno++);
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ",
               "S0_11", vno++, "S022", vno++, "S033",vno++,
               "S0_23", vno++, "S031", vno++, "S012",vno++);
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ",
               "sigma11", vno++, "sigma22", vno++, "sigma33",vno++,
               "sigma23", vno++, "sigma31", vno++, "sigma12",vno++,"sigma_eff",vno++);
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ",
               "S11", vno++, "S22", vno++, "S33",vno++,
               "S23", vno++, "S31", vno++, "S12",vno++);
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] ",
               "E11", vno++, "E22", vno++, "E33",vno++,
               "E23", vno++, "E31", vno++, "E12",vno++,"Eeff",vno++);
  fprintf(out,"%s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d] %s[%d]\n",
               "e11", vno++, "e22", vno++, "e33",vno++,
               "e23", vno++, "e31", vno++, "e12",vno++,"eeff",vno++);  
  return 0;
};

int print_outputs(FILE *out, double t,
                             double J,
                             double dw, double vw, double dX, double vX,
                             double sigma_eff_0, double sigma_eff, double Eeff, double eeff,
                             Matrix(double) *F,  Matrix(double) *sigma0,
                             Matrix(double) *S0, Matrix(double) *sigma,
                             Matrix(double) *S,  Matrix(double) *E,
                             Matrix(double) *e)
{
  fprintf(out,"%e %e ", t, J);
  fprintf(out,"%e %e %e %e %e %e %e %e %e ", 
               Mat_v(*F,1,1),Mat_v(*F,1,2),Mat_v(*F,1,3),
               Mat_v(*F,2,1),Mat_v(*F,2,2),Mat_v(*F,2,3),
               Mat_v(*F,3,1),Mat_v(*F,3,2),Mat_v(*F,3,3));
  fprintf(out,"%e %e %e %e ", dw,vw,dX,vX);
  fprintf(out,"%e %e %e %e %e %e %e ",
               Mat_v(*sigma0,1,1),Mat_v(*sigma0,2,2),Mat_v(*sigma0,3,3),
               Mat_v(*sigma0,2,3),Mat_v(*sigma0,3,1),Mat_v(*sigma0,1,2),sigma_eff_0);
  fprintf(out,"%e %e %e %e %e %e ",
               Mat_v(*S0,1,1),Mat_v(*S0,2,2),Mat_v(*S0,3,3),
               Mat_v(*S0,2,3),Mat_v(*S0,3,1),Mat_v(*S0,1,2));
  fprintf(out,"%e %e %e %e %e %e %e ",
               Mat_v(*sigma,1,1),Mat_v(*sigma,2,2),Mat_v(*sigma,3,3),
               Mat_v(*sigma,2,3),Mat_v(*sigma,3,1),Mat_v(*sigma,1,2),sigma_eff);
  fprintf(out,"%e %e %e %e %e %e ",
               Mat_v(*S,1,1),Mat_v(*S,2,2),Mat_v(*S,3,3),
               Mat_v(*S,2,3),Mat_v(*S,3,1),Mat_v(*S,1,2));
  fprintf(out,"%e %e %e %e %e %e %e ",
               Mat_v(*E,1,1),Mat_v(*E,2,2),Mat_v(*E,3,3),
               Mat_v(*E,2,3),Mat_v(*E,3,1),Mat_v(*E,1,2),Eeff);
  fprintf(out,"%e %e %e %e %e %e %e\n",
               Mat_v(*e,1,1),Mat_v(*e,2,2),Mat_v(*e,3,3),
               Mat_v(*e,2,3),Mat_v(*e,3,1),Mat_v(*e,1,2),eeff);  
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
  int construct_ET = 0; // computing elasticity tensor
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, construct_ET);
    
  // -----------------------
  // construct GL(3) tensors
  // -----------------------
  enum {F,Fnp1,C,E,sigma0,S0,sigma,e,FI,FITE,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_init(double, F2[a],DIM_3,DIM_3,0.0);
    Matrix_eye(F2[a],DIM_3);
  }    

  Matrix(double) S;
  S.m_row = S.m_col = 3; S.m_pdata = elast.S;  // no need free

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
    
    err += (sim->F_of_t)(&(F2[Fnp1]), t);
    
    // compute strains
    Matrix_AxB(F2[C],1.0,0.0,F2[Fnp1],1,F2[Fnp1],0);
    Matrix_AeqB(F2[E],0.5,F2[C]);
    Mat_v(F2[E],1,1) -= 0.5; 
    Mat_v(F2[E],2,2) -= 0.5;
    Mat_v(F2[E],3,3) -= 0.5;
    double Eeff = 0.0;
    Matrix_ddot(F2[E],F2[E],Eeff);
    Eeff = sqrt(Eeff);
    Matrix_inv(F2[Fnp1],F2[FI]);
    Matrix_AxB(F2[FITE], 1.0,0.0,F2[FI],1,F2[E],0);
    Matrix_AxB(F2[e],  1.0,0.0,F2[FITE], 0,F2[FI],1);
    double eeff = 0.0;
    Matrix_ddot(F2[e],F2[e],Eeff);
    eeff = sqrt(eeff);

        
    //compute undamaged values
    double sigma_eff_0 = 0.0;
    
    err += elast.update_elasticity(&elast,F2[Fnp1].m_pdata, construct_ET);
    err += elast.compute_Cauchy_eff(&elast, &sigma_eff_0, F2[Fnp1].m_pdata);
    err += elast.compute_Cauchy(&elast, F2[sigma0].m_pdata, F2[Fnp1].m_pdata);
    Matrix_AeqB(F2[S0],1.0,S);
        
    //compute damaged values
    double sigma_eff = 0.0;
    Matrix_det(F2[Fnp1], Jnp1);

    err += continuum_damage_integration_alg(&mat_d,&elast,
                                            &w,&X,&H,&is_it_damaged,
                                            wn,Xn,dt,F2[Fnp1].m_pdata);
                                                                                 
    err += update_damaged_elasticity(&mat_d,&elast,w,is_it_damaged,H,
                                     dt,F2[Fnp1].m_pdata, construct_ET);
                                            
    err += elast.compute_Cauchy_eff(&elast,&sigma_eff,F2[Fnp1].m_pdata);
    err += elast.compute_Cauchy(&elast,F2[sigma].m_pdata,F2[Fnp1].m_pdata);
    
    err += print_outputs(out,t,Jnp1,w,0.0,X,0.0,
                                  sigma_eff_0,sigma_eff,Eeff,eeff,
                                  &F2[Fnp1], &F2[sigma0],
                                  &F2[S0],   &F2[sigma],
                                  &S,        &F2[E],
                                  &F2[e]);
    wn = w;
    Xn = X;
    is_it_damaged = 0;
  }
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);  

  free(F2);  
    
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
                                             mat->beta_dev,
                                             mat->alpha_vol,
                                             mat->beta_vol);                                             

  //err += print_material_property_elasticity(&mat_e); // just printing
  
  // --------------------------------------------------
  // construct and initialize elasticity handle
  // : provides computing manner of elasticity response
  // --------------------------------------------------
  int construct_ET = 1; // computing elasticity tensor
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, construct_ET);
    
  // -----------------------
  // construct GL(3) tensors
  // -----------------------
  enum {F,Fnp1,C,E,sigma0,S0,sigma,e,FI,FITE,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_init(double, F2[a],DIM_3,DIM_3,0.0);
    Matrix_eye(F2[a],DIM_3);
  }    

  Matrix(double) S;
  S.m_row = S.m_col = 3; S.m_pdata = elast.S;  // no need free

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
    
    err += (sim->F_of_t)(&(F2[Fnp1]), t);
    
    // compute strains
    Matrix_AxB(F2[C],1.0,0.0,F2[Fnp1],1,F2[Fnp1],0);
    Matrix_AeqB(F2[E],0.5,F2[C]);
    Mat_v(F2[E],1,1) -= 0.5; 
    Mat_v(F2[E],2,2) -= 0.5;
    Mat_v(F2[E],3,3) -= 0.5;
    double Eeff = 0.0;
    Matrix_ddot(F2[E],F2[E],Eeff);
    Eeff = sqrt(Eeff);
    Matrix_inv(F2[Fnp1],F2[FI]);
    Matrix_AxB(F2[FITE], 1.0,0.0,F2[FI],1,F2[E],0);
    Matrix_AxB(F2[e],  1.0,0.0,F2[FITE], 0,F2[FI],1);
    double eeff = 0.0;
    Matrix_ddot(F2[e],F2[e],Eeff);
    eeff = sqrt(eeff);

        
    //compute undamaged values
    double sigma_eff_0 = 0.0;
    
    err += elast.update_elasticity(&elast,F2[Fnp1].m_pdata, construct_ET);
    err += elast.compute_Cauchy_eff(&elast, &sigma_eff_0, F2[Fnp1].m_pdata);
    err += elast.compute_Cauchy(&elast, F2[sigma0].m_pdata, F2[Fnp1].m_pdata);
    Matrix_AeqB(F2[S0],1.0,S);
        
    //compute damaged values
    double sigma_eff = 0.0;
    Matrix_det(F2[Fnp1], Jnp1);
    
    err += continuum_damage_split_integration_alg(&mat_d,&elast,
                                     &dw,&vw,&dX,&vX,&dH,&vH,
                                     &is_it_damaged_d,&is_it_damaged_v,                                     
                                     dwn,vwn,dXn,vXn,dt,F2[Fnp1].m_pdata);    
    
    err += update_damaged_elasticity_split(&mat_d,&elast,dw,vw,dH,vH,
                                          is_it_damaged_d,is_it_damaged_v,
                                          dt, F2[Fnp1].m_pdata, construct_ET);
                                          
    err += elast.compute_Cauchy_eff(&elast,&sigma_eff,F2[Fnp1].m_pdata);
    err += elast.compute_Cauchy(&elast,F2[sigma].m_pdata,F2[Fnp1].m_pdata);
    
    err += print_outputs(out,t,Jnp1,dw,vw,dX,vX,
                                  sigma_eff_0,sigma_eff,Eeff,eeff,
                                  &F2[Fnp1], &F2[sigma0],
                                  &F2[S0],   &F2[sigma],
                                  &S,        &F2[E],
                                  &F2[e]);          
    dwn = dw;
    vwn = vw;
    dXn = dX;
    vXn = vX;
    
    is_it_damaged_d = 0;
    is_it_damaged_v = 0;
  }
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);  

  free(F2);  
    
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
  mat.alpha_dev = 1.0;
  mat.beta_dev  = 0.0;
  mat.alpha_vol = 0.0;
  mat.beta_vol  = 1.0;

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

  mat.beta_dev = 0.0;
  sprintf(sim.file_out, "expansion_beta_dev_0p0.txt");
  err += test_split_damage_model(&mat, &sim);  

  mat.beta_dev = 0.3;
  sprintf(sim.file_out, "expansion_beta_dev_0p3.txt");
  err += test_split_damage_model(&mat, &sim);  

  mat.beta_dev = 0.7;
  sprintf(sim.file_out, "expansion_beta_dev_0p7.txt");
  err += test_split_damage_model(&mat, &sim);

  mat.mu = 0.1;
  mat.alpha_vol = 0.0;
  mat.beta_dev = 0.0;  
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
  mat.beta_dev = 0.0;  
  sprintf(sim.file_out, "tension_alpha_vol_0p2.txt");
  err += test_split_damage_model(&mat, &sim);  

  mat.alpha_vol = 0.0;
  mat.beta_dev = 0.2;  
  sprintf(sim.file_out, "tension_beta_dev_0p2.txt");
  err += test_split_damage_model(&mat, &sim);
  
  mat.alpha_vol = 0.2;
  mat.beta_dev = 0.2;  
  sprintf(sim.file_out, "tension_alpha_vol_0p2_beta_dev_0p2.txt");
  err += test_split_damage_model(&mat, &sim);
 
  mat.alpha_vol = 1.0;
  mat.beta_dev = 1.0;  
  sprintf(sim.file_out, "tension_alpha_vol_1p0_beta_dev_1p0.txt");
  err += test_split_damage_model(&mat, &sim);

  sprintf(sim.file_out, "tension_original_damage.txt");
  err += test_damage_model(&mat, &sim);
  
  sim.F_of_t = F_of_t_tension_compression;
  mat.alpha_vol = 0.2;
  mat.beta_dev = 0.2;  
  sprintf(sim.file_out, "tension_compression_beta_dev_0p2_beta_vol_1p0.txt");
  err += test_split_damage_model(&mat, &sim);
  
  mat.beta_vol = -1.0;
  mat.beta_dev = 0.0;  
  sprintf(sim.file_out, "tension_compression_beta_dev_0p0_beta_vol_m1p0.txt");
  err += test_split_damage_model(&mat, &sim);
  
  mat.beta_vol = -1.0;
  mat.beta_dev = -0.2;  
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

