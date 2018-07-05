#include <stdio.h>
#include <math.h>

/// user define total deformation gradients
int F_of_t(double *F,
           double t)
{
  double strainrate = -0.005;
  double   I[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};                  
  double Sh2[9] = {0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};  
   
  double lambdac = (t > 40) ? strainrate*40.0: t*strainrate;
  double lambdas = (t > 40) ? (40.0-t)*strainrate : 0.0;
                   
    for(int ia=0; ia<9; ia++)
      F[ia] = I[ia] + lambdac*I[ia] + lambdas*Sh2[ia];    
  
  return 0;
}           
