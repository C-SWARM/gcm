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

void test_damage_model(void)
{
  int err = 0;
  
  Matrix(double) F,Fnp1,C;
  Matrix_construct_redim(double, F, F_ROW_MAX, F_COLUMN_MAX);
  Matrix_construct_redim(double, Fnp1, 3,3);
  Matrix_construct_redim(double, C, 3,3);
  
  err += read_deformation_gradient(&F);
  
  MATERIAL_ELASTICITY mat_e;
  MATERIAL_CONTINUUM_DAMAGE mat_d;
 
  double p1 = 8.0;
  double p2 = 2.5;
  double Yin = 0.15;
  double mu = 100.0;
 
  set_properties_using_E_and_nu(&mat_e,9.95e+3,0.34);
  set_damage_parameters(&mat_d, p1, p2, Yin, mu, 0.98);
  print_material_property_elasticity(&mat_e);
  
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);
  
  CONTINUUM_DAMAGE_SPLIT damage;

  err += initialize_continuum_damage_split(&damage,&elast,&mat_d,0.0);  
  
  Matrix(double) S, sigma;
  Matrix_construct_redim(double,sigma,3,3); 
  S.m_row = S.m_col = 3; S.m_pdata = elast.S;  // no need free
  
  double dt = 0.001;
  double Jnp1;
  for(int a = 0; a<F_ROW_MAX; a++)
  {
    double t = dt*a;
    
    for(int b = 0; b<9; b++)
      Fnp1.m_pdata[b] = Mat_v(F,a+1,b+2);

    //compute undamaged values
    double sigma_eff_0 = 0.0;
    
    elast.update_elasticity(&elast,Fnp1.m_pdata, 0);
    elast.compute_Cauchy_eff(&elast,&sigma_eff_0,Fnp1.m_pdata);
    elast.compute_Cauchy(&elast,sigma.m_pdata,Fnp1.m_pdata);
    double S11_0 = elast.S[0];
    double sigma_11_0 = sigma.m_pdata[0];
    
    //compute damaged values
    double sigma_eff = 0.0;
    Matrix_det(Fnp1, Jnp1);
    err += continuum_damage_split_integration_alg(&damage, dt, Fnp1.m_pdata);
    update_damaged_elasticity_split(&damage, dt, Fnp1.m_pdata, 0);
    elast.compute_Cauchy_eff(&elast,&sigma_eff,Fnp1.m_pdata);
    elast.compute_Cauchy(&elast,sigma.m_pdata,Fnp1.m_pdata);
    
    Matrix_AxB(C,1.0,0.0,Fnp1,1,Fnp1,0);    
    printf("%e %e %e %e %e %e %e %e %e %e %e %e\n", t, (a+1.0)*0.001, Jnp1, damage.wh, damage.wu, damage.X, 
                                                 sigma_eff_0,       sigma_11_0,      S11_0, 
                                                 sigma_eff,   sigma.m_pdata[0], elast.S[0]);

    err += update_damage_split_time_steps(&damage);

  }       
   
  destruct_elasticity(&elast);
  Matrix_cleanup(F);
  Matrix_cleanup(Fnp1);
  Matrix_cleanup(C);
  Matrix_cleanup(sigma);
}


int main(int argc,char *argv[])
{
  int err = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  test_damage_model();    
  
  gettimeofday(&end, NULL);
  double diff = (double)(end.tv_usec - start.tv_usec)/1000000.0 
  + (double)(end.tv_sec - start.tv_sec);
  printf ("Total time: %.4lf s\n", diff);
  
  return err;
}
