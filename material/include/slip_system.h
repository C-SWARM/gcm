#ifndef H__H__SLIP_SYSTEM__H__H
#define H__H__SLIP_SYSTEM__H__H

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

enum{SLIP_SYSTEM_FCC,SLIP_SYSTEM_BCC,SLIP_SYSTEM_HCP};
enum{EULER_ANGLE_XYZ,EULER_ANGLE_BUNGE};

struct SLIP_SYSTEM
{
  double *p_sys;
  int N_SYS;
  int unit_cell;
  char name[1024];
  int    ort_option[3];
  char   ort_file_in[2048];
  double ort_angles[3];
};

typedef struct SLIP_SYSTEM SLIP_SYSTEM;

int construct_slip_system(SLIP_SYSTEM *slip, int type);
int destruct_slip_system(SLIP_SYSTEM *slip);

/// compute rotation matrix using Euler angles
///
/// \param[out] R_in computed rotation matrxi
/// \param[out] A1_in rotation matrix of angle1
/// \param[out] A2_in rotation matrix of angle2
/// \param[out] A3_in rotation matrix of angle3
/// \param[in] angle1 1st Euler angle
/// \param[in] angle2 2nd Euler angle
/// \param[in] angle3 3rd Euler angle
/// \param[in] type Euler angle type, 0 (default): R = Az*Ay*Az
///                                   1          : Bunge R = Az*Ax*Az 
/// \return non-zero on interal error
int rotation_matrix_of_Euler_angles(double *R_in, 
                                    double *A1_in,
                                    double *A2_in,
                                    double *A3_in,
                                    const double angle1, 
                                    const double angle2, 
                                    const double angle3,
                                    const int type = EULER_ANGLE_XYZ);
                                        
/// compute rotation matrices using array of Euler angles
///
/// \param[out] R_out computed rotation matrices
/// \param[in] angles array of Euler angles angles = [phi_0, theta_0, psi_0
///                                                   phi_1, theta_1, psi_1
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
