#include <stdio.h>
#include <math.h>

/// user define total deformation gradients
int F_of_t(double *F,
           double t)
{
  double d = 0.001;
  double I11[9] = {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double Iii[9] = {0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};  
                      
  for(int ia=0; ia<9; ia++)
    F[ia] = I11[ia]*(1.0 + t*d) + Iii[ia]*(1.0 - 0.5*t*d);
  
  return 0;
}           
