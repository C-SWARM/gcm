#include "constitutive_model.h"
#include "slip_system.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int construct_slip_system(SLIP_SYSTEM *slip, int type)
{
  int err = 0;
  
  double FCC[12*3*3] = {
   0.408248, 0.408248, 0.408248,-0.408248,-0.408248,-0.408248, 0.      , 0.      , 0.      ,
   0.      , 0.      , 0.      ,-0.408248,-0.408248,-0.408248, 0.408248, 0.408248, 0.408248,
   0.408248, 0.408248, 0.408248, 0.      , 0.      , 0.      ,-0.408248,-0.408248,-0.408248,
  -0.408248, 0.408248, 0.408248,-0.408248, 0.408248, 0.408248, 0.      , 0.      , 0.      ,
  -0.408248, 0.408248, 0.408248, 0.      , 0.      , 0.      ,-0.408248, 0.408248, 0.408248,
   0.      , 0.      , 0.      , 0.408248,-0.408248,-0.408248,-0.408248, 0.408248, 0.408248,
   0.408248,-0.408248, 0.408248, 0.408248,-0.408248, 0.408248, 0.      , 0.      , 0.      ,
   0.      , 0.      , 0.      , 0.408248,-0.408248, 0.408248, 0.408248,-0.408248, 0.408248,
  -0.408248, 0.408248,-0.408248, 0.      , 0.      , 0.      , 0.408248,-0.408248, 0.408248,
   0.      , 0.      , 0.      , 0.408248, 0.408248,-0.408248, 0.408248, 0.408248,-0.408248,
   0.408248, 0.408248,-0.408248, 0.      , 0.      , 0.      , 0.408248, 0.408248,-0.408248,
   0.408248, 0.408248,-0.408248,-0.408248,-0.408248, 0.408248, 0.      , 0.      , 0.};
  
  slip->unit_cell = type; 
  switch(type)
  {
    case SLIP_SYSTEM_FCC:
      slip->N_SYS = N_SYS;   //sets to constant declared in slip_system.h
      slip->p_sys = (double *) malloc(sizeof(double)*12*3*3);
      sprintf(slip->name, "FCC");      
      for(int a=0; a<12*3*3; a++)
        slip->p_sys[a] = FCC[a];
      break;        
    case SLIP_SYSTEM_BCC:
    case SLIP_SYSTEM_HCP:
    default:
      break;
  }
  slip->ort_option[0] = -1;
  return err;
}


int destruct_slip_system(SLIP_SYSTEM *slip)
{
  int err = 0;
  if(slip->p_sys!=NULL)
    free(slip->p_sys);
  return err;
}

/// compute rotation matrix using Euler angles
/// R = Az(psi)*Ay(theta)*Ax(phi)
///
/// \param[out] R_in computed rotation matrxi
/// \param[out] Ax_in rotation matrix of phi
/// \param[out] Ay_in rotation matrix of theta
/// \param[out] Az_in rotation matrix of psi
/// \param[in] phi 1st Euler angle
/// \param[in] theta 2nd Euler angle
/// \param[in] psi 3rd Euler angle
/// \return non-zero on interal error
int rotation_matrix_of_Euler_angles(double *R_in, 
                                    double *Ax_in,
                                    double *Ay_in,
                                    double *Az_in,
                                    double phi, 
                                    double theta, 
                                    double psi)
{
  int err = 0;
  // handle for matrix objecs based on [3x3] matrix stencil 
  Matrix(double) Ax, Ay, Az, R;
  Ax.m_row = Ax.m_col = DIM_3; Ax.m_pdata = Ax_in;
  Ay.m_row = Ay.m_col = DIM_3; Ay.m_pdata = Ay_in;
  Az.m_row = Az.m_col = DIM_3; Az.m_pdata = Az_in;
   R.m_row =  R.m_col = DIM_3;  R.m_pdata =  R_in;

  // compute rotation matrix of the 1st Euler angle phi
  Mat_v(Ax,1,1) = 1.0;
  Mat_v(Ax,1,2) = Mat_v(Ax,1,3) = Mat_v(Ax,2,1) = Mat_v(Ax,3,1) = 0.0;
  Mat_v(Ax,2,2) =  cos(phi); Mat_v(Ax,2,3) = -sin(phi);
  Mat_v(Ax,3,2) =  sin(phi); Mat_v(Ax,3,3) = cos(phi);    

  // compute rotation matrix of the 2nd Euler angle theta
  Mat_v(Ay,2,2) = 1.0;
  Mat_v(Ay,2,1) = Mat_v(Ay,2,3) = Mat_v(Ay,1,2) = Mat_v(Ay,3,2) = 0.0;
  Mat_v(Ay,1,1) = cos(theta); Mat_v(Ay,1,3) =  sin(theta);
  Mat_v(Ay,3,1) = -sin(theta); Mat_v(Ay,3,3) =  cos(theta);

  // compute rotation matrix of the 3rd Euler angle psi
  Mat_v(Az,3,3) = 1.0;
  Mat_v(Az,3,1) = Mat_v(Az,3,2) = Mat_v(Az,1,3) = Mat_v(Az,2,3) = 0.0;
  Mat_v(Az,1,1) =  cos(psi); Mat_v(Az,1,2) = -sin(psi);
  Mat_v(Az,2,1) =  sin(psi); Mat_v(Az,2,2) = cos(psi);
  
  // R = Az(psi)*Ay(theta)*Ax(phi)            
  Matrix_Tns2_AxBxC(R,1.0,0.0,Az,Ay,Ax);
   
  return err;
}

/// compute rotation matrices using array of Euler angles
///
/// \param[out] R_out computed rotation matrices
/// \param[in] angles array of Euler angles angles = [phi_0, theta_0, psi_0
///                                                    phi_1, theta_1, psi_1
///                                                              :
///                                                                   ]
/// \param[in] ortno number of orientation angles
/// \return non-zero on interal error 
int set_crystal_orientations(double *R_out, double *angles, int ortno)
{
  int err = 0;
 
  double Ax[DIM_3x3], Az[DIM_3x3], Ay[DIM_3x3]; 
    
  for(int a = 0; a<ortno; a++)
  {
    double phi   = angles[a*DIM_3 + 0];
    double theta = angles[a*DIM_3 + 1];
    double psi   = angles[a*DIM_3 + 2];
    
    err += rotation_matrix_of_Euler_angles(R_out + DIM_3x3*a, 
                                           Ax, Ay, Az, phi, theta, psi);    
  }
  return err;
}


int rotate_crystal_orientation(SLIP_SYSTEM* slip_out, double *R_in, SLIP_SYSTEM* slip_in)
{
  int err = 0;
  Matrix(double) R, RT, Pa_in, Pa_out;

   Pa_in.m_row =  Pa_in.m_col = DIM_3;
  Pa_out.m_row = Pa_out.m_col = DIM_3;
       R.m_row =      R.m_col = DIM_3; R.m_pdata = R_in;           
  
  Matrix_construct_redim(double, RT, DIM_3,DIM_3);
  Matrix_AeqBT(RT,1.0,R);
  
  for(int a=0; a<slip_in->N_SYS; a++)
  {
     Pa_in.m_pdata =  (slip_in->p_sys) + a*DIM_3x3;
    Pa_out.m_pdata = (slip_out->p_sys) + a*DIM_3x3;
    Matrix_Tns2_AxBxC(Pa_out,1.0,0.0,R,Pa_in,RT);    
  }
    
  Matrix_cleanup(RT);
  return err;    
}

int generate_random_crystal_orientation(double* angles, int angle_no)
{
  int err = 0;
  double pi = 3.141592653589793;
  double psi_max = pi;

  time_t t;  
  srand((unsigned) time(&t));

  for(int a=0; a<angle_no; a++)
  {
    double n1 = (double)rand() / ((double)RAND_MAX + 1);
    double n2 = (double)rand() / ((double)RAND_MAX + 1);
    double n3 = (double)rand() / ((double)RAND_MAX + 1);
    
    double phi   = 2.0*pi*n1;
    double theta = asin(n2);
    double psi   = psi_max*(2.0*n3-1.0);
    
    angles[a*DIM_3+0] = phi;      
    angles[a*DIM_3+1] = theta;
    angles[a*DIM_3+2] = psi;            
  }
  return err;
} 
