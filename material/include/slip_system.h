#ifndef H__H__SLIP_SYSTEM__H__H
#define H__H__SLIP_SYSTEM__H__H

enum{SLIP_SYSTEM_FCC,SLIP_SYSTEM_BCC,SLIP_SYSTEM_HCP};

struct SLIP_SYSTEM
{
  double *p_sys;
  int N_SYS;
  char name[1024];
};

typedef struct SLIP_SYSTEM SLIP_SYSTEM;

int construct_slip_system(SLIP_SYSTEM *slip, int type);
int destruct_slip_system(SLIP_SYSTEM *slip);

int set_crystal_orientations(double *R_out, double *angles, int ortno);
int rotate_crystal_orientation(SLIP_SYSTEM* slip_out, double *R, SLIP_SYSTEM* slip_in);
int generate_random_crystal_orientation(double* angles, int angle_no);
#endif