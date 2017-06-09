#include "hardening.h"
#include "material_properties.h"
#include <math.h>

double compute_g_np1(double gs_np1, double gn, double dt, double gamma_dot_np1, 
                     MATERIAL_CRYSTAL_PLASTICITY *mat_p)
{
  double A = ((gs_np1-mat_p->g0)*gn + dt*mat_p->G0*gs_np1*gamma_dot_np1);
  double B = ((gs_np1-mat_p->g0) + dt*mat_p->G0*gamma_dot_np1);
  return A/B; // g_np1
} 

double compute_gs_np1(double gamma_dot_np1, MATERIAL_CRYSTAL_PLASTICITY *mat_p)
{
   return mat_p->gs_0*pow(fabs(gamma_dot_np1/mat_p->gamma_dot_s), mat_p->w); // gs_np1
}
