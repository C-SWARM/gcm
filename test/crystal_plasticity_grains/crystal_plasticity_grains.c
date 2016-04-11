#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "crystal_plasticity_integration.h"
#include "flowlaw.h"
//#include <sys/time.h>
#include "mpi.h"

#define VNO 18

typedef int (*DEFORMATION_GRADIENT) (Matrix(double) *Fn, Matrix(double) *F, Matrix(double) *L, double t, double dt);

int F_of_t_tension(Matrix(double) *Fn, 
                       Matrix(double) *F, 
                       Matrix(double) *L_in,
                       double t, double dt)
{
  int err = 0;
  Matrix(double) L;
  Matrix_construct_init(double,L,DIM_3,DIM_3,0.0);

  double d = 1.0;
  L.m_pdata[0] = d;
  L.m_pdata[4] = L.m_pdata[8] = -d*0.5;
  
  Matrix_AeqB(*L_in, 1.0, L);
  Fnp1_Implicit(F->m_pdata, Fn->m_pdata, L.m_pdata, dt);  
    
  Matrix_cleanup(L);
  return err;
}; 

int F_of_t_relaxation(Matrix(double) *Fn, 
                       Matrix(double) *F, 
                       Matrix(double) *L_in,
                       double t, double dt)
{
  int err = 0;
  Matrix(double) L;
  Matrix_construct_init(double,L,DIM_3,DIM_3,0.0);

  if(t<=0.5)
  {  
    double d = 1.0;
    L.m_pdata[0] = -d;
    L.m_pdata[4] = L.m_pdata[8] = d*0.5;
  }
  
  Matrix_AeqB(*L_in, 1.0, L);
  Fnp1_Implicit(F->m_pdata, Fn->m_pdata, L.m_pdata, dt);  
    
  Matrix_cleanup(L);
  return err;
}; 

int F_of_t_cyclic(Matrix(double) *Fn, 
                  Matrix(double) *F, 
                  Matrix(double) *L_in,
                  double t, double dt)
{
  int err = 0;
  Matrix(double) L;
  Matrix_construct_init(double,L,DIM_3,DIM_3,0.0);

  if(t<=0.5)
  {  
    double d = 1.0;
    L.m_pdata[0] = d;
    L.m_pdata[4] = L.m_pdata[8] = -d*0.5;
  }
  if((0.5 < t) && (t <= 0.9))
  {
    double d = 1.0;
    L.m_pdata[0] = -d;
    L.m_pdata[4] = L.m_pdata[8] = d*0.5;
  }
  if(0.9<t)
  {
    double d = 1.0;
    L.m_pdata[0] = d;
    L.m_pdata[4] = L.m_pdata[8] = -d*0.5;
  }  
  
  Matrix_AeqB(*L_in, 1.0, L);
  Fnp1_Implicit(F->m_pdata, Fn->m_pdata, L.m_pdata, dt);  
    
  Matrix_cleanup(L);
  return err;
};                        
typedef struct
{
  double lame1;
  double lame2;
  double E;
  double nu;

  double gamma_dot_0;
  double gamma_dot_s;
  double m;
  double g0;
  double G0;
  double gs_0;
  double w;
  int slip_system;
} MAT_PROP; 

typedef struct {
  double dt;
  int stepno;
  char file_out[1024];
  DEFORMATION_GRADIENT F_of_t;
} SIM_PARAMS;

void test_crystal_plasticity_single_crystal(MAT_PROP *mat_in, SIM_PARAMS *sim,
                                            double *result_out, double *angle)
{
  
  int max_itr_stag      = 100;
  int max_itr_hardening = 5;
  int max_itr_M         = 100;
  double tol_hardening  = 1.0e-6;
  double tol_M          = 1.0e-6;
  double computer_zero  = 1.0e-15;

  // create material properties: Elasticity
  MATERIAL_ELASTICITY mat_e;
  if(mat_in->E<0)  
    set_properties_using_Lame_constants(&mat_e,mat_in->lame1,mat_in->lame2);
  else
    set_properties_using_E_and_nu(&mat_e,mat_in->E,mat_in->nu);
  
  mat_e.devPotFlag = 2;
  mat_e.volPotFlag = 99;     
  //print_material_property_elasticity(&mat_e); // <= this is optional
  
  // create slip system : 0 for FCC and rotate  
  SLIP_SYSTEM slip_in, slip;
  construct_slip_system(&slip_in,mat_in->slip_system);   
  construct_slip_system(&slip,   mat_in->slip_system);
  
  double R[DIM_3x3];
  set_crystal_orientations(R, angle, 1);
  rotate_crystal_orientation(&slip, R, &slip_in);
    
  destruct_slip_system(&slip_in); 
  
  // create material properties: Plasticity
  MATERIAL_CRYSTAL_PLASTICITY mat_p;
  set_properties_crystal_plasticity(&mat_p,&slip,
                                     mat_in->gamma_dot_0,
                                     mat_in->gamma_dot_s, 
                                     mat_in->m,
                                     mat_in->g0,
                                     mat_in->G0,
                                     mat_in->gs_0,
                                     mat_in->w);                                     
  //print_material_property_crystal_plasticity(&mat_p);  // <= this is optional 

  // create material plasticity: it needs material properties for elasticity and plasticity
  MATERIAL_CONSTITUTIVE_MODEL mat;
  set_properties_constitutive_model(&mat,&mat_e,&mat_p);
  
  // create solver info: criteria for numerical iterations
  CRYSTAL_PLASTICITY_SOLVER_INFO solver_info;
  set_crystal_plasticity_solver_info(&solver_info,max_itr_stag,
                                                  max_itr_hardening,
                                                  max_itr_M,
                                                  tol_hardening,
                                                  tol_M,
                                                  computer_zero);
  //solver_info.max_subdivision = 128;                                                  
  //print_crystal_plasticity_solver_info(&solver_info); // <= this is optional
  
  // create elasticity object for integration
  // this creates memory for stress and elasticity tensor s.t. requires destructor
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);  

  // set variables for integration
  enum {M,MI,pFn,pFnp1,pFnp1_I,eFnp1,Fn,Fnp1,L,D,sigma,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_init(double, F2[a],DIM_3,DIM_3,0.0);
    Matrix_eye(F2[a],DIM_3);
  } 
  
  double g_n,g_np1;
  g_n = g_np1 = mat_p.g0;  
  
  // start integration  
  Matrix(double) result;
  result.m_row = sim->stepno;
  result.m_col = VNO;
  result.m_pdata = result_out;  

  for(int a = 1; a<=sim->stepno; a++)
  {
    double lambda = 0.0;
    double t = a*(sim->dt);
    
    // compute total deformation gradient using velocity gradient
    sim->F_of_t(F2+Fn,F2+Fnp1,F2+L,t,sim->dt); 
    
    staggered_Newton_Rapson(F2[pFnp1].m_pdata,F2[M].m_pdata, &g_np1, &lambda, 
                            F2[pFn].m_pdata, F2[Fn].m_pdata,F2[Fnp1].m_pdata, 
                            g_n, sim->dt, &mat, &elast, &solver_info);

    Matrix_AeqB(F2[pFn],1.0,F2[pFnp1]);
    Matrix_AeqB(F2[Fn],1.0,F2[Fnp1]);  
    Matrix_inv(F2[pFnp1],F2[pFnp1_I]);
    Matrix_AxB(F2[eFnp1],1.0,0.0,F2[Fnp1],0,F2[pFnp1_I],0);    
    
    g_n = g_np1;
    
    // print result at time t    
    elast.update_elasticity(&elast,F2[eFnp1].m_pdata,0);
    elast.compute_Cauchy(&elast,F2[sigma].m_pdata,F2[eFnp1].m_pdata);

    Matrix_symmetric(F2[L],F2[D]);

    for(int ij=1; ij<=DIM_3x3; ij++)
    {
      Mat_v(result, a, ij)           += F2[sigma].m_pdata[ij-1];
      Mat_v(result, a, ij + DIM_3x3) +=     F2[D].m_pdata[ij-1];
    }        
  }    
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);  

  free(F2);    
  destruct_elasticity(&elast);
  destruct_slip_system(&slip);
}


int test_crystal_plasticity_grains(MAT_PROP *mat,
                                    SIM_PARAMS *sim, 
                                    int myrank, 
                                    int nproc, 
                                    int n_grain)
{
  int err = 0;
  Matrix(double) angle, result, G_result;
  Matrix_construct_init(double, angle,   n_grain,     DIM_3, 0.0);
  Matrix_construct_init(double, result,  sim->stepno,   VNO, 0.0);
  Matrix_construct_init(double, G_result,sim->stepno,   VNO, 0.0);    
  
  // generate random angles -->
  if(myrank==0)
    generate_random_crystal_orientation(angle.m_pdata, n_grain);

  MPI_Bcast(angle.m_pdata,n_grain*DIM_3, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
  // <-- generate random angles
    
  if(nproc>=n_grain && myrank<n_grain)
  {  
    double *temp = (angle.m_pdata) + myrank*DIM_3;
    test_crystal_plasticity_single_crystal(mat,sim,result.m_pdata,temp);
    printf("grain[%3d] is computed.\n", myrank);              
  } 

  if(nproc<n_grain)
  {
    int dNG = (int) n_grain/nproc;
    for(int a=0; a<dNG+1; a++)
    {
      if(a*nproc + myrank>=n_grain)
        break;

      double *temp = (angle.m_pdata) + (a*nproc + myrank)*DIM_3;                
      test_crystal_plasticity_single_crystal(mat,sim,result.m_pdata,temp);
      printf("grain[%3d] is computed.\n", a*nproc + myrank);                          
    }
  }
  
  MPI_Reduce(result.m_pdata,G_result.m_pdata,(sim->stepno)*VNO,
             MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  if(myrank==0)
  {
    Matrix(double) T, D;
    T.m_row = T.m_col = DIM_3;
    D.m_row = D.m_col = DIM_3;

    char fname[1024];
    sprintf(fname, "%s_%d.txt", sim->file_out,n_grain); 
    FILE *fp = fopen(fname, "w");
    
    double e_eff = 0.0;
    double e_11 = 0.0;
    for(int a=0; a<(sim->stepno)*VNO; a++)
      G_result.m_pdata[a] = G_result.m_pdata[a]/n_grain;
      
    for(int a=0; a<sim->stepno; a++)
    {
      T.m_pdata = G_result.m_pdata + a*VNO;
      D.m_pdata = G_result.m_pdata + (a*VNO + DIM_3x3);
      
      double T_eff = (Mat_v(T,1,1) - Mat_v(T,2,2))*(Mat_v(T,1,1) - Mat_v(T,2,2)) + 
                     (Mat_v(T,2,2) - Mat_v(T,3,3))*(Mat_v(T,2,2) - Mat_v(T,3,3)) + 
                     (Mat_v(T,3,3) - Mat_v(T,1,1))*(Mat_v(T,3,3) - Mat_v(T,1,1)) +
                     6.0*(Mat_v(T,1,2)*Mat_v(T,1,2)+
                          Mat_v(T,2,3)*Mat_v(T,2,3)+
                          Mat_v(T,3,1)*Mat_v(T,3,1));

      T_eff = sqrt(T_eff/2.0);
        
      double DD = 0.0;
      Matrix_ddot(D,D,DD);
      e_eff += sqrt(2.0/3.0*DD)*sim->dt;
      e_11 += D.m_pdata[0]*sim->dt;
      fprintf(fp, "%e %e %e %e %e\n",(a+1)*sim->dt,e_eff,T_eff, e_11, T.m_pdata[0]);
    }
    fclose(fp);
  }  
  		      
  Matrix_cleanup(angle);
  Matrix_cleanup(result);
  Matrix_cleanup(G_result);
  return 0;    
}

int main(int argc,char *argv[])
{
  int err = 0;
  
  int myrank = 0;
  int nproc = 0;
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank (mpi_comm,&myrank);
  MPI_Comm_size (mpi_comm,&nproc);
    
  int n_grain = 1;
  if(argc>=2)
  {  
    sscanf(argv[1], "%d", &n_grain);
    if(n_grain<=0)
      n_grain = 1; 
  }
  
  if(myrank==0)
    printf("Number of grains : %d\n", n_grain);

  double total_time = 0.0;
  total_time -= MPI_Wtime();
  
  MAT_PROP mat_1998_Dawson = {
    .lame1       = 75.6e+3,
    .lame2       = 26.1e+3,
    .E           = -1, //will be computed later
    .nu          = -1, //
    .gamma_dot_0 = 1.0,
    .gamma_dot_s = 50.0e+9,
    .m           = 0.05,  
    .g0          = 210.0,
    .G0          = 200.0,
    .gs_0        = 330.0,
    .w           = 0.005,
  };

  MAT_PROP mat_1999_Zacharia = {
    .lame1       = 54.4e+3,
    .lame2       = 25.3e+3,
    .E           = -1, //will be computed later
    .nu          = -1, //
    .gamma_dot_0 = 1.0,
    .gamma_dot_s = 50.0e+9,
    .m           = 0.05,  
    .g0          = 27.17,
    .G0          = 58.41,
    .gs_0        = 61.8,
    .w           = 0.005,
  };  

  SIM_PARAMS sim = {
    .dt = 0.0001,
    .stepno = 20000,
    .F_of_t = F_of_t_tension,
  };
  
  sprintf(sim.file_out, "tension");
  err += test_crystal_plasticity_grains(&mat_1999_Zacharia, &sim, myrank, nproc, n_grain);    

  sprintf(sim.file_out, "relaxation");
  sim.dt = 0.001;
  sim.F_of_t = F_of_t_relaxation;
  sim.stepno = 1500;
  err += test_crystal_plasticity_grains(&mat_1998_Dawson, &sim, myrank, nproc, n_grain);

  sprintf(sim.file_out, "cyclic");
  sim.dt = 0.0001;
  sim.F_of_t = F_of_t_cyclic;
  sim.stepno = 16000;
  err += test_crystal_plasticity_grains(&mat_1998_Dawson, &sim, myrank, nproc, n_grain);

  total_time += MPI_Wtime();
  if(myrank==0) printf ("Total time: %.4lf s\n", total_time);
  
  MPI_Finalize();
  return err;
}
