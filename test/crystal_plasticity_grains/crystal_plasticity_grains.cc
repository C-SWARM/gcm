#include "constitutive_model.h"
#include "material_properties.h"
#include "hyperelasticity.h"
#include "crystal_plasticity_integration.h"
#include "flowlaw.h"
#include "mpi.h"

#define VNO 18
#define MAX_GRAIN 1024

typedef struct {
  double dt;             /// time step size
  int    stepno;         /// number of time steps
  char   file_out[1024]; /// output file name
  int    intg_type;      /// 
  double t1;
  double t2;
  Matrix(double) L1;
  Matrix(double) L2;
} SIM_PARAMS;

/// compute total deformation gradient by integrating velocity gradient
///
/// \param[in]  *Fn  deformation gradient at tn
/// \param[out] *F   deformation gradient at t
/// \param[out] *L   velocity gradient
/// \param[in]   t   current time
/// \param[in]  *sim simulation parameters are defined e.g dt and time step size and velocity gradient 
/// \return non-zero on internal error
int F_of_t(Matrix(double) *Fn,
           Matrix(double) *F,
           Matrix(double) *L,
           double t,            
           const SIM_PARAMS *sim)
{
  int err = 0;
  
  Matrix_init(*L, 0.0);
  switch(sim->intg_type)
  {
    case 1:
      Matrix_AeqB(*L, 1.0, sim->L1);
      break;
    case 2:
      if(t<=sim->t1)
        Matrix_AeqB(*L, 1.0, sim->L1);

      if((sim->t1 < t) && (t <= sim->t2))
        Matrix_AeqB(*L, 1.0, sim->L2);

      if(sim->t2<t)
        Matrix_AeqB(*L, 1.0, sim->L1);
      
      break;  
    case 3:
      if(t<=sim->t1)
        Matrix_AeqB(*L, 1.0, sim->L1);
      
      break;
    default:
      Matrix_AeqB(*L, 1.0, sim->L1);
      break;
  }      
  
  Fnp1_Implicit(F->m_pdata, Fn->m_pdata, L->m_pdata, sim->dt);  
  return err;
};

typedef struct
{
  double lame1;       /// 1st leme constant
  double lame2;       /// 2nd leme constant
  double E;           /// Young's modulus
  double nu;          /// Poisson's ratio
  double gamma_dot_0; /// initial shearing rate
  double gamma_dot_s; /// shearing rate saturation
  double m;           /// rate sensitivity
  double g0;          /// initial hardening
  double G0;          /// initial shear strength
  double gs_0;        /// material parameter for hardening saturation
  double w;           /// material parameter for hardening saturation
  int slip_system;    /// 0: FCC, 1: BCC, or 2: HCP
} MAT_PROP; 

/// perform integration algorithm for a single grain
///
/// \param[in]  *mat_in     material parameters
/// \param[in]  *sim        simulation parameters
/// \param[out] *result_out array of results
/// \param[in]  *angle      list of Euler angle
/// \param[in]  *grain_id   grain ID
int test_crystal_plasticity_single_crystal(const MAT_PROP *mat_in, 
                                           const SIM_PARAMS *sim,
                                           double *result_out, 
                                           double *angle,
                                           const int grain_id)
{
  int err = 0;
  int max_itr_stag      = 100;
  int max_itr_hardening = 5;
  int max_itr_M         = 100;
  double tol_hardening  = 1.0e-6;
  double tol_M          = 1.0e-6;
  double computer_zero  = 1.0e-15;

  // create material properties: Elasticity
  MATERIAL_ELASTICITY mat_e;
  if(mat_in->E<0)  
    err += set_properties_using_Lame_constants(&mat_e,mat_in->lame1,mat_in->lame2);
  else
    err += set_properties_using_E_and_nu(&mat_e,mat_in->E,mat_in->nu);
  
  mat_e.devPotFlag = 1;
  mat_e.volPotFlag = 2;     
  // print_material_property_elasticity(&mat_e); // <= this is optional
  // create slip system : 0 for FCC and rotate  
  SLIP_SYSTEM slip_in, slip;
  err += construct_slip_system(&slip_in,mat_in->slip_system);   
  err += construct_slip_system(&slip,   mat_in->slip_system);
  
  double R[DIM_3x3];
  err += set_crystal_orientations(R, angle, 1);
  err += rotate_crystal_orientation(&slip, R, &slip_in);
    
  err += destruct_slip_system(&slip_in); 
  
  // create material properties: Plasticity
  MATERIAL_CRYSTAL_PLASTICITY mat_p;
  err += set_properties_crystal_plasticity(&mat_p,&slip,
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
  err += set_properties_constitutive_model(&mat,&mat_e,&mat_p);
  
  // create solver info: criteria for numerical iterations
  CRYSTAL_PLASTICITY_SOLVER_INFO solver_info;
  err += set_crystal_plasticity_solver_info(&solver_info,max_itr_stag,
                                             max_itr_hardening,
                                             max_itr_M,
                                             tol_hardening,
                                             tol_M,
                                             computer_zero);
  solver_info.max_subdivision = 128;                                                  
  //print_crystal_plasticity_solver_info(&solver_info); // <= this is optional
  
  // create elasticity object for integration
  // this creates memory for stress and elasticity tensor s.t. requires destructor
  ELASTICITY elast;
  err += construct_elasticity(&elast, &mat_e, 1);  

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
    F_of_t(F2+Fn,F2+Fnp1,F2+L,t,sim); 
    
    err += staggered_Newton_Rapson(F2[pFnp1].m_pdata,F2[M].m_pdata, &g_np1, &lambda, 
                                   F2[pFn].m_pdata, F2[Fn].m_pdata,F2[Fnp1].m_pdata, 
                                   g_n, sim->dt, &mat, &elast, &solver_info);
    
    Matrix_AeqB(F2[pFn],1.0,F2[pFnp1]);
    Matrix_AeqB(F2[Fn],1.0,F2[Fnp1]);  
    Matrix_inv(F2[pFnp1],F2[pFnp1_I]);
    Matrix_AxB(F2[eFnp1],1.0,0.0,F2[Fnp1],0,F2[pFnp1_I],0);    
    
    g_n = g_np1;
    
    // print result at time t    
    err += elast.update_elasticity(&elast,F2[eFnp1].m_pdata,0);
    err += elast.compute_Cauchy(&elast,F2[sigma].m_pdata,F2[eFnp1].m_pdata);

    Matrix_symmetric(F2[L],F2[D]);

    for(int ij=1; ij<=DIM_3x3; ij++)
    {
      Mat_v(result, a, ij)           += F2[sigma].m_pdata[ij-1];
      Mat_v(result, a, ij + DIM_3x3) +=     F2[D].m_pdata[ij-1];
    }        
  }
  
  char out[1024];
  sprintf(out, "%s_Fs_%d.txt", sim->file_out, grain_id);
  FILE *fp = fopen(out, "w");
  fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
              F2[pFnp1].m_pdata[0], F2[pFnp1].m_pdata[1], F2[pFnp1].m_pdata[2],
              F2[pFnp1].m_pdata[3], F2[pFnp1].m_pdata[4], F2[pFnp1].m_pdata[5],
              F2[pFnp1].m_pdata[6], F2[pFnp1].m_pdata[7], F2[pFnp1].m_pdata[8],
              F2[Fnp1].m_pdata[0], F2[Fnp1].m_pdata[1], F2[Fnp1].m_pdata[2],
              F2[Fnp1].m_pdata[3], F2[Fnp1].m_pdata[4], F2[Fnp1].m_pdata[5],
              F2[Fnp1].m_pdata[6], F2[Fnp1].m_pdata[7], F2[Fnp1].m_pdata[8]);
  fclose(fp);
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);  

  free(F2);    
  err += destruct_elasticity(&elast);
  err += destruct_slip_system(&slip);
  return err;
}

/// read orientations
/// 
/// Read Euler angles from file
///
/// \param[out] *angles  list of angles read from file
/// \param[out] *n_grain number of grains read from file
/// \param[in]  *fn_in   filebase name to be read, filename = [fn_in]_[myrank].in
/// \param[in]  myrank  partion ID of the orientation file
int read_orientations(Matrix(double) *angles, int *n_grain, char *fn_in, int myrank)
{
  int err = 0;
  char fn[1024], line[1024];
  sprintf(fn, "%s_%d.in", fn_in, myrank);
  FILE *fp = fopen(fn, "r");
  if(fp==NULL)
  { 
    printf("fail to read [%s]\n", fn);
    err++;
    return err;
  }
  
  int cnt = 0;
  while(fgets(line, 1024, fp)!=NULL)
  {
    if(line[0]=='#')
	    continue;
        
    int e, ip;    
    double x1, x2, x3;        
    sscanf(line, "%d %d %lf %lf %lf", &e, &ip, &x1, &x2, &x3);
    int id = cnt + *n_grain;
    Mat_v(*angles, id+1, 1) = x1;    
    Mat_v(*angles, id+1, 2) = x2;        
    Mat_v(*angles, id+1, 3) = x3;
    cnt++;
  } 

  *n_grain += cnt;
  
  fclose(fp);
  return err;
}

/// Spawn single integration algorithm as may as the number of grains
///
/// First, rank 0 process read orientation angles and number of grains 
/// and broadcat the read values to other processes. 
/// Each process with the grain angles, takes uniformly distributed grains 
/// and perform single grain integrations at a time. 
/// When every process is done with their jobs, results are MPI_reduced to write outputs. 
///
/// \param[in] *mat        material properties
/// \param[in] *simulation parameters
/// \param[in] myrank      MPI rank
/// \param[in] nproc       number of processes 
/// \param[in] *fn_in      filebase name of orientations
/// \return non-zero on internal errror 
int test_crystal_plasticity_grains(MAT_PROP *mat,
                                   SIM_PARAMS *sim, 
                                   int myrank, 
                                   int nproc, 
                                   char *fn_in,
                                   int NP)
{
  int err = 0;

  Matrix(double) angle, result, G_result;
  
  // read orientation angle -->
  int cnt_file_read = 0;
  int n_grain = 0;

  if(myrank==0)
  {
    Matrix(double) temp_angle;
    Matrix_construct_init(double, temp_angle, MAX_GRAIN, DIM_3, 0.0);
    for(int ia=0; ia<NP; ia++)
    {
      err += read_orientations(&temp_angle, &n_grain, fn_in, ia);
      if(err>0)
        break;
    }
    
    if(err==0)
    {
      Matrix_construct_init(double, angle, n_grain, DIM_3, 0.0);
      for(int ia=0; ia<n_grain*DIM_3; ia++)
        angle.m_pdata[ia] = temp_angle.m_pdata[ia];
    }
    Matrix_cleanup(temp_angle);
  }
  
  MPI_Bcast(&err,1, MPI_INT, 0, MPI_COMM_WORLD);

  if(err>0)
  {
    printf("Error on reading crystal orientations\n");
    return err;
  }
  
  MPI_Bcast(&n_grain,1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if(myrank>0)
    Matrix_construct_init(double, angle, n_grain, DIM_3, 0.0);

  MPI_Bcast(angle.m_pdata,n_grain*DIM_3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if(myrank==0)
  {
    FILE *fp_angle = fopen("angles.txt", "w");
    for(int ia=1; ia<=n_grain; ia++)
    {
      for(int ib=1; ib<=DIM_3; ib++)
        fprintf(fp_angle, "%e ", Mat_v(angle, ia, ib));
      fprintf(fp_angle, "\n");
    }
    fclose(fp_angle);
  }
  // <-- read orientation angle

  Matrix_construct_init(double, result,  sim->stepno,   VNO, 0.0);
  Matrix_construct_init(double, G_result,sim->stepno,   VNO, 0.0);
      
  if(nproc>=n_grain && myrank<n_grain)
  {  
    double *temp = (angle.m_pdata) + myrank*DIM_3;
    err += test_crystal_plasticity_single_crystal(mat,sim,result.m_pdata,temp,myrank);
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
      err += test_crystal_plasticity_single_crystal(mat,sim,result.m_pdata,temp,a*nproc + myrank);
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
  
  char fn_sim[1024];
  char fn_ort[1024];  
  int NP = -1;
  
  MAT_PROP mat;
  SIM_PARAMS sim;
  double L1[9];
  double L2[9];

  sim.L1.m_row = sim.L1.m_col = DIM_3; sim.L1.m_pdata = L1;
  sim.L2.m_row = sim.L2.m_col = DIM_3; sim.L2.m_pdata = L2;    

  if(argc<4)
  {
    if(myrank==0)
    {  
      printf("Usage: mpirun -np # ./crystal_plasticity_grains [FILE_SIM] [FILE_ORT] [np]\n");
      printf("\t[FILE_SIM]\t: file path with simulation parameters\n");
      printf("\t[FILE_ORT]\t: file path for orientation, e.g CO/co\n");
      printf("\t          \t: CO contains co_0.in\n");
      printf("\t          \t:             co_1.in\n");
      printf("\t          \t:                :   \n");
      printf("\t          \t:        co_(np-1).in\n");
      printf("\t[np]\t\t: number of partitions of the orientation files\n");              
    }
    exit(-1);
  }
  else
  {
    // read commend line arguments
    sscanf(argv[1], "%s", fn_sim);
    sscanf(argv[2], "%s", fn_ort);    
    sscanf(argv[3], "%d", &NP);
    if(NP<=0)
    {
      if(myrank==0)
        printf("Number of partitions is smaller than 1 (NP=%d).\nExit\n", NP);
      
      return -1;
    }
    
    // read material and simulation parameters
    // if fail to read the file, exit
    FILE *fp_sim = NULL;
    fp_sim = fopen(fn_sim, "r");
    
    if(fp_sim==NULL)
    { 
      if(myrank==0)
        printf("fail to read [%s]\n", fn_sim);

      return -1;      
    }
    
    char line[1024]; 
    // read material parameters 
    while(fgets(line, 1024, fp_sim)!=NULL)
    {
      if(line[0]!='#')
        break;
    }    
    sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
                 &(mat.lame1), &(mat.lame2), &(mat.E), &(mat.nu), 
                 &(mat.gamma_dot_0), &(mat.gamma_dot_s), &(mat.m), 
                 &(mat.g0), &(mat.G0), &(mat.gs_0), &(mat.w), &(mat.slip_system));

    // read analysis name
    while(fgets(line, 1024, fp_sim)!=NULL)
    {
      if(line[0]!='#')
        break;
    }
    sscanf(line, "%s", sim.file_out);
    
    // read time steps
    while(fgets(line, 1024, fp_sim)!=NULL)
    {
      if(line[0]!='#')
        break;
    }
    sscanf(line, "%lf %d", &(sim.dt), &(sim.stepno));
    
    // read integration option
    while(fgets(line, 1024, fp_sim)!=NULL)
    {
      if(line[0]!='#')
        break;
    }    
    sscanf(line, "%d", &(sim.intg_type));

    // read 1st velocity gradient
    while(fgets(line, 1024, fp_sim)!=NULL)
    {
      if(line[0]!='#')
        break;
    }

    sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", L1+0,L1+1,L1+2,
                                                            L1+3,L1+4,L1+5,
                                                            L1+6,L1+7,L1+8,
                                                            &(sim.t1));
    if(sim.intg_type==2)
    {
      // read 2nd velocity gradient
      while(fgets(line, 1024, fp_sim)!=NULL)
      {
        if(line[0]!='#')
          break;
      }      
                                                                
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", L2+0,L2+1,L2+2,
                                                              L2+3,L2+4,L2+5,
                                                              L2+6,L2+7,L2+8, 
                                                              &(sim.t2));
    }
    fclose(fp_sim);
  }

  double total_time = 0.0;
  total_time -= MPI_Wtime();

  err += test_crystal_plasticity_grains(&mat, &sim, myrank, nproc, fn_ort, NP);

  total_time += MPI_Wtime();
  if(myrank==0) printf ("Total time: %.4lf s\n", total_time);
  
  MPI_Finalize();
  return err;
}
