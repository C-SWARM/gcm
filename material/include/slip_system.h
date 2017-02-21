#ifndef H__H__SLIP_SYSTEM__H__H
#define H__H__SLIP_SYSTEM__H__H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

const int N_SYS = 12;  //helpful for declaring ttl tensors of size slip->N_SYS

enum{SLIP_SYSTEM_FCC,SLIP_SYSTEM_BCC,SLIP_SYSTEM_HCP};

struct SLIP_SYSTEM
{
  double *p_sys;
  int N_SYS;
  int unit_cell;
  char name[1024];
  int    ort_option[2];
  char   ort_file_in[2048];
  double ort_angles[3];
};

typedef struct SLIP_SYSTEM SLIP_SYSTEM;

int construct_slip_system(SLIP_SYSTEM *slip, int type);
int destruct_slip_system(SLIP_SYSTEM *slip);

/// compute rotation matrix using Euler angles
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
                                    double psi);

/// compute rotation matrices using array of Euler angles
///
/// \param[out] R_out computed rotation matrices
/// \param[in] angles array of Euler angles angles = [phi_0, theta_0, psi_0
///                                                    phi_1, theta_1, psi_1
///                                                              :
///                                                                   ]
/// \param[in] ortno number of orientation angles
/// \return non-zero on interal error 
int set_crystal_orientations(double *R_out, double *angles, int ortno);

int rotate_crystal_orientation(SLIP_SYSTEM* slip_out, double *R, SLIP_SYSTEM* slip_in);
int generate_random_crystal_orientation(double* angles, int angle_no);


#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
