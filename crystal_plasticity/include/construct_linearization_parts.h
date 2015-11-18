#ifndef H__H__CONSTRUCT_LINEARIZATION_PARTS__H__H
#define H__H__CONSTRUCT_LINEARIZATION_PARTS__H__H

int compute_Kuu_a(double *Kuu_a_out, double *A_in,  double *S_in, 
                  double *C_in,      double *Pa_in, double *L_in, double drdtaus_a);
                                    
int compute_Kuu_b(double *Kuu_b_out, double *MI_in);

#endif                  