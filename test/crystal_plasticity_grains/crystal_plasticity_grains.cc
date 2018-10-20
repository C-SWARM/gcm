#include "constitutive_model_handle.h"
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
  Tensor<2> L1;
  Tensor<2> L2;
} SIM_PARAMS;

/// compute total deformation gradient by integrating velocity gradient
///
/// \param[in]  *Fn  deformation gradient at tn
/// \param[out] *F   deformation gradient at t
/// \param[out] *L   velocity gradient
/// \param[in]   t   current time
/// \param[in]  *sim simulation parameters are defined e.g dt and time step size and velocity gradient 
/// \return non-zero on internal error
template <class T1, class T2, class T3>
int F_of_t(T1 &Fn,
           T2 &F,
           T3 &L,
           double t,            
           const SIM_PARAMS *sim)
{
  int err = 0;
  
  for(int ia=0; ia<DIM_3x3; ia++)
    L.data[ia] = 0.0;
    
  switch(sim->intg_type)
  {
    case 1:
      L(i,j) = sim->L1(i,j);
      break;
    case 2:
      if(t<=sim->t1)
        L(i,j) = sim->L1(i,j);
      if((sim->t1 < t) && (t <= sim->t2))
        L(i,j) = sim->L2(i,j);

      if(sim->t2<t)
        L(i,j) = sim->L1(i,j);
      
      break;  
    case 3:
      if(t<=sim->t1)
        L(i,j) = sim->L1(i,j);
      
      break;
    default:
      L(i,j) = sim->L1(i,j);
      break;
  }      
  
  Fnp1_Implicit(F.data, Fn.data, L.data, sim->dt);  
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
  GcmSolverInfo solver_info;
  int max_sub_cycling = 128;
  set_gcm_solver_info(&solver_info,max_itr_stag,
                      max_itr_hardening,
                      max_itr_M,
                      tol_hardening,
                      tol_M,
                      computer_zero,
                      max_sub_cycling);   

  //print_crystal_plasticity_solver_info(&solver_info); // <= this is optional
  
  // create elasticity object for integration
  // this creates memory for stress and elasticity tensor s.t. requires destructor
  HyperElasticity elast;
  elast.construct_elasticity(&mat_e, true);  

  // set variables for integration
  Tensor<2> pFn,pFnp1,pFnp1_I,eFnp1,Fnm1,Fn,Fnp1,L,D,sigma;

  pFn   = ttl::identity(i,j);
  pFnp1 = ttl::identity(i,j);
  eFnp1 = ttl::identity(i,j);
  Fnm1  = ttl::identity(i,j);
  Fn    = ttl::identity(i,j);
  Fnp1  = ttl::identity(i,j);
    
  double lambda = 0.0;
  double g_np1;
  
  double hF[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};  
  
  GcmCpIntegrator integrator;
  
  integrator.mat         = &mat;
  integrator.elasticity  = &elast;
  integrator.solver_info = &solver_info; 
    
  integrator.set_tensors(Fnp1.data,
                         Fn.data,
                         Fnm1.data,
                         pFnp1.data,
                         pFn.data,
                         hF,hF);
                         
  integrator.gnp1   = &g_np1;
  integrator.lambda = &lambda;
  integrator.gn     = g_np1 = mat_p.g0;
  
  for(int a = 1; a<=sim->stepno; a++)
  {
    double t = a*(sim->dt);
    lambda = 0.0;
    
    // compute total deformation gradient using velocity gradient
    F_of_t(Fn,Fnp1,L,t,sim); 
    
    err += integrator.run_integration_algorithm(sim->dt, sim->dt);
       
    pFn  = pFnp1(i,j);
    Fnm1 = Fn(i,j);
    Fn   = Fnp1(i,j);
    err += inv(pFnp1,pFnp1_I);
    eFnp1 = Fnp1(i,k)*pFnp1_I(k,j);
         
    integrator.gn = g_np1;
    
    // print result at time t    
    err += elast.update_elasticity(eFnp1.data,false);
    elast.compute_Cauchy(sigma.data,eFnp1.data);

    D(i,j) = 0.5*(L(i,j) + L(j,i));

    for(int ia=0; ia<DIM_3x3; ia++)
    {
      result_out[VNO*(a-1) + ia]           += sigma.data[ia];
      result_out[VNO*(a-1) + DIM_3x3 + ia] +=     D.data[ia];
    }        
  }
  
  char out[2048];
  sprintf(out, "%s_Fs_%d.txt", sim->file_out, grain_id);
  FILE *fp = fopen(out, "w");
  fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
              pFnp1.data[0], pFnp1.data[1], pFnp1.data[2],
              pFnp1.data[3], pFnp1.data[4], pFnp1.data[5],
              pFnp1.data[6], pFnp1.data[7], pFnp1.data[8],
               Fnp1.data[0],  Fnp1.data[1],  Fnp1.data[2],
               Fnp1.data[3],  Fnp1.data[4],  Fnp1.data[5],
               Fnp1.data[6],  Fnp1.data[7],  Fnp1.data[8]);
  fclose(fp);
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
int read_orientations(double *angles, int *n_grain, char *fn_in, int myrank)
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
    if(id>MAX_GRAIN)
    {  
      printf("Read grains more than maximum number of grains (MAX_GRAIN = %d)\nExit.", MAX_GRAIN);
      exit(-1); 
    }
    angles[id*DIM_3+0] = x1;
    angles[id*DIM_3+1] = x2;
    angles[id*DIM_3+2] = x3;
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

  double *angle    = NULL;
  double *result   = NULL;
  double *G_result = NULL;
  
  // read orientation angle -->
  int n_grain = 0;

  if(myrank==0)
  {
    double *temp_angle = new double[MAX_GRAIN*DIM_3];
 
    for(int ia=0; ia<NP; ia++)
    {
      err += read_orientations(temp_angle, &n_grain, fn_in, ia);
      if(err>0)
        break;
    }
    
    if(err==0)
    {
      angle = new double[n_grain*DIM_3];
      for(int ia=0; ia<n_grain*DIM_3; ia++)
        angle[ia] = temp_angle[ia];
    }
    delete temp_angle;
  }
  
  MPI_Bcast(&err,1, MPI_INT, 0, MPI_COMM_WORLD);

  if(err>0)
  {
    printf("Error on reading crystal orientations\n");
    return err;
  }
    
  MPI_Bcast(&n_grain,1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if(myrank>0)
    angle = new double[n_grain*DIM_3];

  MPI_Bcast(angle,n_grain*DIM_3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if(myrank==0)
  {
    FILE *fp_angle = fopen("angles.txt", "w");
    for(int ia=0; ia<n_grain; ia++)
    {
      for(int ib=0; ib<DIM_3; ib++)
        fprintf(fp_angle, "%e ", angle[ia*DIM_3+ib]);
      fprintf(fp_angle, "\n");
    }
    fclose(fp_angle);
  }
  // <-- read orientation angle

    result = new double[sim->stepno*VNO]();
  G_result = new double[sim->stepno*VNO]();
      
  if(nproc>=n_grain && myrank<n_grain)
  {  
    double *temp = angle + myrank*DIM_3;
    err += test_crystal_plasticity_single_crystal(mat,sim,result,temp,myrank);
    printf("grain[%3d] is computed.\n", myrank);              
  } 

  if(nproc<n_grain)
  {
    int dNG = (int) n_grain/nproc;
    for(int a=0; a<dNG+1; a++)
    {
      if(a*nproc + myrank>=n_grain)
        break;

      double *temp = angle + (a*nproc + myrank)*DIM_3;                
      err += test_crystal_plasticity_single_crystal(mat,sim,result,temp,a*nproc + myrank);
      printf("grain[%3d] is computed.\n", a*nproc + myrank);                          
    }
  }
  
  MPI_Reduce(result,G_result,(sim->stepno)*VNO,
             MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  if(myrank==0)
  {
    char fname[2048];
    sprintf(fname, "%s_%d.txt", sim->file_out,n_grain); 
    FILE *fp = fopen(fname, "w");

    double e_eff = 0.0;
    double e_11 = 0.0;
    for(int a=0; a<(sim->stepno)*VNO; a++)
      G_result[a] = G_result[a]/n_grain;
      
    for(int a=0; a<sim->stepno; a++)
    {
      TensorA<2> T(G_result + a*VNO);
      TensorA<2> D(G_result + a*VNO + DIM_3x3);
      
      double T_eff = (T[0][0] - T[1][1])*(T[0][0] - T[1][1]) + 
                     (T[1][1] - T[2][2])*(T[1][1] - T[2][2]) + 
                     (T[2][2] - T[0][0])*(T[2][2] - T[0][0]) +
                     6.0*(T[0][1]*T[0][1]+
                          T[1][2]*T[1][2]+
                          T[2][0]*T[2][0]);

      T_eff = sqrt(T_eff/2.0);
        
      double DD = D(i,j)*D(i,j);
      e_eff += sqrt(2.0/3.0*DD)*sim->dt;
      e_11 += D[0][0]*sim->dt;
      fprintf(fp, "%e %e %e %e %e\n",(a+1)*sim->dt,e_eff,T_eff, e_11, T[0][0]);
    }
    fclose(fp);
  } 
  
  delete angle;
  delete result;
  delete G_result;
  		      
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
  double *L1 = sim.L1.data;
  double *L2 = sim.L2.data;

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
    mat.slip_system = 0;
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
