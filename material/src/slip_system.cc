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
      slip->N_SYS = 12;
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
  slip->ort_option[2] = EULER_ANGLE_XYZ;
  return err;
}


int destruct_slip_system(SLIP_SYSTEM *slip)
{
  int err = 0;
  if(slip->p_sys!=NULL)
    free(slip->p_sys);
  return err;
}

/// compute rotation matrix using Euler angles (Roe rotation)
/// R = Az(psi)*Ay(theta)*Ax(phi)
///
/// \param[out] R_in computed rotation matrxi
/// \param[out] Ax_in rotation matrix of phi
/// \param[out] Ay_in rotation matrix of theta
/// \param[out] Az_in rotation matrix of psi
/// \param[in] phi 1st Euler angle
/// \param[in] theta 2nd Euler angle
/// \param[in] psi 3rd Euler angle
void rotation_matrix_of_Euler_angles_Roe(double *R_in, 
                                         double *Ax_in,
                                         double *Ay_in,
                                         double *Az_in,
                                         const double phi, 
                                         const double theta, 
                                         const double psi)
{
  // handle for matrix objecs based on [3x3] matrix stencil 
  TensorA<2> R(R_in), Ax(Ax_in), Ay(Ay_in), Az(Az_in);
  
  // compute rotation matrix of the 1st Euler angle phi
  Ax[0][0] = 1.0;
  Ax[0][1] = Ax[0][2] = Ax[1][0] = Ax[2][0] = 0.0;
  Ax[1][1] =  cos(phi); Ax[1][2] = -sin(phi);
  Ax[2][1] =  sin(phi); Ax[2][2] =  cos(phi);    

  // compute rotation matrix of the 2nd Euler angle theta
  Ay[1][1] = 1.0;
  Ay[1][0] = Ay[1][2] = Ay[0][1] = Ay[2][1] = 0.0;
  Ay[0][0] =  cos(theta); Ay[0][2] = sin(theta);
  Ay[2][0] = -sin(theta); Ay[2][2] = cos(theta);

  // compute rotation matrix of the 3rd Euler angle psi
  Az[2][2] = 1.0;
  Az[2][0] = Az[2][1] = Az[0][2] = Az[1][2] = 0.0;
  Az[0][0] =  cos(psi); Az[0][1] = -sin(psi);
  Az[1][0] =  sin(psi); Az[1][1] =  cos(psi);
  
  R = Az(i,j)*Ay(j,k)*Ax(k,l);
}

/// compute rotation matrix using Euler angles (Bunge)
/// R = Az(psi2)*Ax(phi)*Az(psi1)
///
/// \param[out] R_in computed rotation matrxi
/// \param[out] Az1_in rotation matrix of psi1
/// \param[out] Ax_in  rotation matrix of phi
/// \param[out] Az2_in rotation matrix of psi2
/// \param[in] psi1 1st Euler angle
/// \param[in] phi  2nd Euler angle
/// \param[in] psi2 3rd Euler angle
/// \return non-zero on interal error
void rotation_matrix_of_Euler_angles_Bunge(double *R_in, 
                                           double *Az1_in,
                                           double *Ax_in,
                                           double *Az2_in,
                                           const double psi1, 
                                           const double phi, 
                                           const double psi2)
{
  // handle for matrix objecs based on [3x3] matrix stencil 
  TensorA<2> R(R_in), Az1(Az1_in), Ax(Ax_in), Az2(Az2_in);

  // compute rotation matrix of the 1st Euler angle psi_1
  Az1[2][2] = 1.0;
  Az1[2][0] = Az1[2][1] = Az1[0][2] = Az1[1][2] = 0.0;
  Az1[0][0] =  cos(psi1);  Az1[0][1] = sin(psi1);
  Az1[1][0] = -sin(psi1);  Az1[1][1] = cos(psi1);

  // compute rotation matrix of the 2nd Euler angle phi
  Ax[0][0] = 1.0;
  Ax[0][1] = Ax[0][2] = Ax[1][0] = Ax[2][0] = 0.0;
  Ax[1][1] =  cos(phi);  Ax[1][2] = sin(phi);
  Ax[2][1] = -sin(phi);  Ax[2][2] = cos(phi);    

  // compute rotation matrix of the 3rd Euler angle psi_2
  Az2[2][2] = 1.0;
  Az2[2][0] = Az2[2][1] = Az2[0][2] = Az2[1][2] = 0.0;
  Az2[0][0] =  cos(psi2);  Az2[0][1] = sin(psi2);
  Az2[1][0] = -sin(psi2);  Az2[1][1] = cos(psi2);
  
  R = Az2(i,j)*Ax(j,k)*Az1(k,l);
}

/// compute rotation matrix using Euler angles
///
/// \param[out] R computed rotation matrxi
/// \param[out] A1 rotation matrix of angle1
/// \param[out] A2 rotation matrix of angle2
/// \param[out] A3 rotation matrix of angle3
/// \param[in] angle1 1st Euler angle
/// \param[in] angle2 2nd Euler angle
/// \param[in] angle3 3rd Euler angle
/// \param[in] type Euler angle type, 0 (default): R = Az*Ay*Az
///                                   1          : Bunge R = Az*Ax*Az 
/// \return non-zero on interal error
int rotation_matrix_of_Euler_angles(double *R, 
                                    double *A1,
                                    double *A2,
                                    double *A3,
                                    const double angle1, 
                                    const double angle2, 
                                    const double angle3,
                                    const int type)
{
  switch(type)
  {
    case EULER_ANGLE_XYZ:
      rotation_matrix_of_Euler_angles_Roe(R, A1, A2, A3, angle1, angle2, angle3);
      break;
    case EULER_ANGLE_BUNGE:
      rotation_matrix_of_Euler_angles_Bunge(R, A1, A2, A3, angle1, angle2, angle3);
      break;
    default:
      rotation_matrix_of_Euler_angles_Roe(R, A1, A2, A3, angle1, angle2, angle3);
  }
  return 0;
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
  TensorA<2> R(R_in);
    
  for(int a=0; a<slip_in->N_SYS; a++)
  {
    TensorA<2> Pa((slip_in->p_sys) + a*DIM_3x3);
    TensorA<2> Pa_out((slip_out->p_sys) + a*DIM_3x3);
    Pa_out = R(i,j)*Pa(j,k)*R(l,k);
  }
    
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
