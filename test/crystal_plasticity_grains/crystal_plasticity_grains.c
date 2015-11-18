#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "solve_system.h"
#include "flowlaw.h"
//#include <sys/time.h>
#include "mpi.h"

#define VNO 6

void test_crystal_plasticity_single_crystal(double *L, double *result_out, double *angle,
                                            double dt, double t_step_no)
{
  double lame1 = 75600.0;
  double lame2     = 26100.0;
  double E = 70.0e+3;
  double nu = 0.25;
  
  double gamma_dot_0 = 1.0;
  double gamma_dot_s = 50.0e+9;
  double m           = 0.05;  
  double g0          = 210.0;
  double G0          = 200.0;
  double gs_0        = 330.0;
  double w           = 0.005;
  
  int max_itr_stag      = 100;
  int max_itr_hardening = 5;
  int max_itr_M         = 100;
  double tol_hardening  = 1.0e-6;
  double tol_M          = 1.0e-6;
  double computer_zero  = 1.0e-15;

  // create material properties: Elasticity
  MATERIAL_ELASTICITY mat_e;  
  set_properties_using_E_and_nu(&mat_e,E,nu); 
  // or you can use : set_properties_using_Lame_constants(&mat_e,lame1,lame2);
  //print_material_property_elasticity(&mat_e); // <= this is optional
  
  // create slip system : 0 for FCC and rotate
  // it is needed to setup plasticity
  SLIP_SYSTEM slip_in, slip;
  construct_slip_system(&slip_in,0);   
  construct_slip_system(&slip,0);
  
  double R[DIM_3x3];
  set_crystal_orientations(R, angle, 1);
  rotate_crystal_orientation(&slip, R, &slip_in);
    
  destruct_slip_system(&slip_in); 
  
  
  // create material properties: Plasticity
  MATERIAL_CRYSTAL_PLASTICITY mat_p;
  set_properties_crystal_plasticity(&mat_p,&slip,gamma_dot_0,gamma_dot_s, 
                                     m,g0,G0,gs_0,w);
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
  //print_crystal_plasticity_solver_info(&solver_info); // <= this is optional
  
  // create elasticity object for integration
  // this creates memory for stress and elasticity tensor s.t. requres destructor
  ELASTICITY elast;
  construct_elasticity(&elast, &mat_e, 1);  

  // set variables for integration
  enum {M,MI,pFn,pFnp1,Fn,Fnp1,sigma,PK2dev,sigma_dev,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_init(double, F2[a],DIM_3,DIM_3,0.0);
    Matrix_eye(F2[a],DIM_3);
  } 
  
  double g_n,g_np1;
  g_n = g_np1 = mat_p.g0;  
  
  // start integration  
  Matrix(double) PK2, result;
  PK2.m_row = PK2.m_col = DIM_3; PK2.m_pdata = elast.S;
  result.m_row = t_step_no;
  result.m_col = VNO;
  result.m_pdata = result_out;  
    
  for(int a = 1; a<=t_step_no; a++)
  {
    double lambda = 0.0;
    double t = a*dt;
    
    // compute total deformation gradient using velocity gradient
    Fnp1_Implicit(F2[Fnp1].m_pdata, F2[Fn].m_pdata, L, dt); 
    
    staggered_Newton_Rapson(F2[pFnp1].m_pdata,F2[M].m_pdata, &g_np1, &lambda, 
                            F2[pFn].m_pdata, F2[Fn].m_pdata,F2[Fnp1].m_pdata, 
                            g_n, dt, &mat, &elast, &solver_info);
    Matrix_AeqB(F2[pFn],1.0,F2[pFnp1]);
    Matrix_AeqB(F2[Fn],1.0,F2[Fnp1]);  
    
    g_n = g_np1;
    
    
    // print result at time t
    double trPK2, tr_sigma;
    
    Matrix_trace(PK2,trPK2);
    Matrix_trace(F2[sigma],tr_sigma);
    Matrix_eye(F2[PK2dev], DIM_3);
    Matrix_eye(F2[sigma_dev], DIM_3);
        
    Matrix_AplusB(F2[PK2dev],    1.0, PK2,      -trPK2/3.0, F2[PK2dev]);
    Matrix_AplusB(F2[sigma_dev], 1.0, F2[sigma], -tr_sigma/3.0, F2[sigma_dev]);    
    
    double norm_sigma, norm_PK2;
    Matrix_ddot(F2[PK2dev],F2[PK2dev],norm_PK2);    
    Matrix_ddot(F2[sigma_dev],F2[sigma_dev],norm_sigma);
    
    double sigma_eff=sqrt(3.0/2.0*norm_sigma);
    double PK2_eff = sqrt(3.0/2.0*norm_PK2);    

    Mat_v(result, a, 1) = t;
    Mat_v(result, a, 2) += sigma_eff;
    Mat_v(result, a, 3) += PK2_eff;
    Mat_v(result, a, 4) += g_np1;    
    Mat_v(result, a, 5) += 0.0;
    Mat_v(result, a, 6) += Mat_v(PK2,1,1);       
  }    
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);  

  free(F2);    
  destruct_elasticity(&elast);
  destruct_slip_system(&slip);
}


void test_crystal_plasticity_grains(int myrank, int nproc, int n_grain)
{

  double dt = 0.001;
  int t_step_no = 1000;
          
  Matrix(double) angle, result, G_result, L;
  Matrix_construct_init(double, angle,   n_grain, DIM_3, 0.0);
  Matrix_construct_init(double, result,  t_step_no, VNO, 0.0);
  Matrix_construct_init(double, G_result,t_step_no, VNO, 0.0);    
  Matrix_construct_init(double, L,        DIM_3, DIM_3, 0.0);    
  
  double d = 1.0;
  // set velocity gradient  
  Mat_v(L,1,1) = -d;
  Mat_v(L,2,2) = Mat_v(L,3,3) = d/2;
  
  // generate random angles -->
  if(myrank==0)
    generate_random_crystal_orientation(angle.m_pdata, n_grain);

  MPI_Bcast(angle.m_pdata,n_grain*DIM_3, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
  // <-- generate random angles
    
  if(nproc>=n_grain && myrank<n_grain)
  {  
    double *temp = (angle.m_pdata) + myrank*DIM_3;
    test_crystal_plasticity_single_crystal(L.m_pdata,result.m_pdata,temp,
                                           dt,t_step_no);
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
      test_crystal_plasticity_single_crystal(L.m_pdata,result.m_pdata,temp,
                                             dt, t_step_no);
      printf("grain[%3d] is computed.\n", a*nproc + myrank);                          
    }
  }
  
  MPI_Reduce(result.m_pdata,G_result.m_pdata,t_step_no*VNO,
             MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  if(myrank==0)
  {
    char fname[1024];
    sprintf(fname, "crystal_results_grains_%d.txt", n_grain); 
    FILE *fp = fopen(fname, "w");
    for(int a=1; a<t_step_no; a++)
    {
      fprintf(fp, "%e %e %e %e %e %e\n",Mat_v(result,a,1),         Mat_v(G_result,a,2)/n_grain,
                                        Mat_v(G_result,a,3)/n_grain, Mat_v(G_result,a,4)/n_grain,
                                        Mat_v(G_result,a,5)/n_grain, Mat_v(G_result,a,6)/n_grain);
    }
    fclose(fp);
  }  
  		      
  Matrix_cleanup(angle);
  Matrix_cleanup(result);
  Matrix_cleanup(L);  
  
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
    
  test_crystal_plasticity_grains(myrank, nproc, n_grain);  
  
  total_time += MPI_Wtime();
  if(myrank==0) printf ("Total time: %.4lf s\n", total_time);
  
  MPI_Finalize();
  return err;
}
