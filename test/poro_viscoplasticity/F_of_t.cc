#include <stdio.h>
#include <math.h>

/// user define total deformation gradients
int F_of_t(double *F,
           double t)
{
  double strainrate = -0.005;
  double   I[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};                  
  double Sh2[9] = {0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};  
   
  double lambdac = (t > 42.6) ?strainrate * 42.6*sin(42.6*M_PI/66.0):strainrate*t*sin(t*M_PI/66.0);
  double lambdas = (t > 42.6) ? ( 42.6 - t ) * strainrate : 0.0;
                   
    for(int ia=0; ia<9; ia++)
      F[ia] = I[ia] + lambdac*I[ia] + lambdas*Sh2[ia];    
  
  return 0;
}           