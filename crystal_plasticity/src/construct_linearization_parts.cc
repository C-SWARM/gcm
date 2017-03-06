#include "constitutive_model.h"
#include "construct_linearization_parts.h"

int compute_Kuu_a(double *Kuu_a_out, double *A_in,  double *S_in, 
                  double *C_in,      double *Pa_in, double *L_in, double drdtaus_a)
{
  Matrix(double) Kuu_a, A, S, C, Pa, L;
  
  Kuu_a.m_row = Kuu_a.m_col = DIM_3x3; Kuu_a.m_pdata = Kuu_a_out;
      A.m_row =     A.m_col = DIM_3;       A.m_pdata =  A_in;
      S.m_row =     S.m_col = DIM_3;       S.m_pdata =  S_in;
      C.m_row =     C.m_col = DIM_3;       C.m_pdata =  C_in;
     Pa.m_row =    Pa.m_col = DIM_3;      Pa.m_pdata = Pa_in;
      L.m_row = DIM_3x3x3x3; L.m_col = 1;  L.m_pdata =  L_in;
  
  int err = 0;  
  Matrix_init(Kuu_a, 0.0);

  enum {AA,L_C,Dtau_dM,F2end};
  Matrix(double) *F2 = malloc(F2end*sizeof(Matrix(double)));
  for (int a = 0; a < F2end; a++) {
    Matrix_construct_redim(double, F2[a],DIM_3,DIM_3);    
  } 
  
  Matrix_init(F2[AA],0.0);
  Matrix_AxB(F2[AA],1.0,1.0,Pa,0,S,0);
  Matrix_AxB(F2[AA],1.0,1.0,S,0,Pa,1);

  Matrix_Tns4_dd_Tns2(F2[L_C], L, C);
  Matrix_AxB(F2[AA],1.0,1.0,F2[L_C],0,Pa,0);
  Matrix_init(F2[Dtau_dM],0.0);
  Matrix_AxB(F2[Dtau_dM],1.0,0.0,A,1,F2[AA],0);
  
  for(int a=1; a<=DIM_3; a++)
  {
    for(int b=1; b<=DIM_3; b++)
    {
      int id_a = (a-1)*DIM_3 + b;
      for(int c=1; c<=DIM_3; c++)
      {
        for(int d=1; d<=DIM_3; d++)
        {
          int id_b = (c-1)*DIM_3+d;
          Mat_v(Kuu_a,id_a,id_b) = drdtaus_a*Mat_v(F2[Dtau_dM],a,b)*Mat_v(Pa,c,d); 
        }  
      }
    }
  }  
  
  for(int a = 0; a < F2end; a++)
    Matrix_cleanup(F2[a]);     
  
  free(F2);
   return err;
}

int compute_Kuu_b(double *Kuu_b_out, double *MI_in)
{
  int err = 0;
  Matrix(double) Kuu_b, MI;
  Kuu_b.m_row = Kuu_b.m_col = DIM_3x3; Kuu_b.m_pdata = Kuu_b_out;
     MI.m_row =    MI.m_col = DIM_3;      MI.m_pdata = MI_in;
        
  Matrix_init(Kuu_b, 0.0);
  
  for(int I=1; I<=DIM_3; I++)
  {
    for(int J=1; J<=DIM_3; J++)
    {
      int id_a = (I-1)*DIM_3 + J;
      for(int K=1; K<=DIM_3; K++)
      {
        for(int L=1; L<=DIM_3; L++)
        {
          int id_b = (L-1)*DIM_3 + K;
          Mat_v(Kuu_b,id_a,id_b) = Mat_v(MI,J,I)*Mat_v(MI,K,L) + Mat_v(MI,K,I)*Mat_v(MI,J,L);
        }
      }
    }
  }

 return err;
}
