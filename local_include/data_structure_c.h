#ifndef H__H__datastructure_c__H__H
#define H__H__datastructure_c__H__H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "mkl_cblas.h"
#include "math_help.h"

extern void dgesv( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                double* b, int* ldb, int* info );

/**
 * Define a Matrix structure of type T.
 */
#define Define_Matrix(T)                                                \
typedef struct Matrix_##T                                               \
{                                                                       \
  long m_row, m_col;                                                    \
  long sizeof_T;                                                        \
  T *m_pdata;                                                           \
  T *temp;                                                              \
} Matrix_##T                                                            \

/**
 * Template-class-like declaration of a Matrix structure of type T.
 */
#define Matrix(T) Matrix_##T

/**
 * Indexing functions. !!NOTE!! Indexing starts from 1
 *
 * Provides direct read/write access to data.
 */
#define Vec_v(p, m) (p).m_pdata[(m)-1]
#define Mat_v(p, m, n) (p).m_pdata[((m) - 1)*(p).m_col+(n) - 1]
#define Tns4_v(p, I_,J_,K_,L_) (p).m_pdata[((I_) - 1)*3*3*3 + ((J_) - 1)*3*3 + ((K_) - 1)*3 + (L_) - 1]

/**
 * Set the row/col size of a 4th-order tensor. Does not allocate any
 * memory.
 */
#define Matrix_Tns4_mat_9x9(p) do{                                      \
  (p).m_row = 9;                                                        \
  (p).m_col = 9;                                                        \
} while(0)

/**
 * Convert a Matrix to a (row) Vector
 */
#define Matrix_Mat2Vec(p) do{                                           \
  (p).m_row = (p).m_row*(p).m_col;                                      \
  (p).m_col = 1;                                                        \
} while(0)

/**
 * Convert a (row) Vector to a [m x n] Matrix
 */
#define Matrix_Vec2Mat(p, m, n) do{                                     \
  (p).m_row = (m);                                                      \
  (p).m_col = (n);                                                      \
} while(0)

/**
 * Compute the 4th order identity tensor
 */
#define Matrix_Tns4_eye(p) do{                                          \
  Matrix_init(p,0.0);                                                   \
  for(int I=1; I<=3; I++)                                               \
  {                                                                     \
    for(int J=1; J<=3; J++)                                             \
    {                                                                   \
      for(int K=1; K<=3; K++)                                           \
      {                                                                 \
        for(int L=1; L<=3; L++)                                         \
          Tns4_v(p, I,J,K,L) = 1.0*(I==K)*(J==L);                       \
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while(0)

/**
 * Comptue the 4th order symmetric identity tensor
 */
#define Matrix_Tns4_II(p) do{                                           \
  Matrix_init(p,0.0);                                                   \
  for(int I=1; I<=3; I++)                                               \
  {                                                                     \
    for(int J=1; J<=3; J++)                                             \
    {                                                                   \
      for(int K=1; K<=3; K++)                                           \
      {                                                                 \
        for(int L=1; L<=3; L++)                                         \
          Tns4_v(p, I,J,K,L) = 0.5*((I==K)*(J==L)+(I==L)*(J==K));       \
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while(0)

/**
 * See Matrix_Tns4_eye, transpose k,l
 */
#define Matrix_Tns4_eye_bar(p) do{                                      \
  Matrix_init(p,0.0);                                                   \
  for(int I=1; I<=3; I++)                                               \
  {                                                                     \
    for(int J=1; J<=3; J++)                                             \
    {                                                                   \
      for(int K=1; K<=3; K++)                                           \
      {                                                                 \
        for(int L=1; L<=3; L++)                                         \
          Tns4_v(p, I,J,K,L) = 1.0*(I==L)*(J==K);                       \
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while(0)


#define Matrix_construct(T, p) do {                                     \
  (p).m_pdata = NULL;                                                   \
  (p).temp    = NULL;                                                   \
  (p).m_row = 0;                                                        \
  (p).m_col = 0;                                                        \
  (p).sizeof_T = sizeof(T);                                             \
} while(0) 

#define Matrix_construct_redim(T, p, m, n) do {                         \
  (p).temp    = NULL;                                                   \
  (p).m_row = m;                                                        \
  (p).m_col = n;                                                        \
  (p).sizeof_T = sizeof(T);                                             \
  (p).m_pdata = (T *) malloc((p).sizeof_T*(m)*(n));                     \
} while(0)

#define Matrix_redim(p, m, n) do {                                      \
  Matrix_cleanup(p);                                                    \
  (p).m_pdata =malloc((p).sizeof_T*(m)*(n));                            \
  (p).temp    = NULL;                                                   \
  (p).m_row = m;                                                        \
  (p).m_col = n;                                                        \
} while(0)                                     

#define Matrix_check_null_and_redim(p, m, n) do {                       \
  if((p).m_pdata==NULL || (p).m_row != m || (p).m_col != n)             \
    Matrix_redim(p, m, n);                                              \
} while(0)

#define Matrix_init(p, value) do {                                      \
  long a__A;                                                            \
  for(a__A = 0; a__A < (p).m_row*(p).m_col; a__A++){                    \
      (p).m_pdata[a__A] =  value;                                       \
  }                                                                     \
} while(0) 

#define Matrix_construct_init(T, p, m, n, value) do {                   \
  (p).temp    = NULL;                                                   \
  (p).m_row = m;                                                        \
  (p).m_col = n;                                                        \
  (p).sizeof_T = sizeof(T);                                             \
  (p).m_pdata = malloc((p).sizeof_T*(m)*(n));                           \
/*  memset(p.m_pdata,value,p.sizeof_T*m*n);   */                        \
  Matrix_init(p, value);                                                \
} while(0)
                                                

#define Matrix_init_w_array(p, m, n, q) do {                            \
  Matrix_check_null_and_redim(p, m, n);                                 \
  long a__A, b__B;                                                      \
  for(a__A = 1; a__A <= (p).m_row; a__A++){                             \
    for(b__B = 1; b__B <= (p).m_col; b__B++){                           \
      Mat_v(p, a__A, b__B) =  (q)[(a__A-1)*(n) + (b__B-1)];             \
    }                                                                   \
  }                                                                     \
} while(0)

/* A = delta_ij */
#define Matrix_eye(A, m) do {                                           \
  long a__A;                                                            \
  Matrix_check_null_and_redim(A, m, m);                                 \
  Matrix_init(A, 0.0);                                                  \
  for(a__A = 1; a__A <= m; a__A++)                                      \
    Mat_v(A, a__A, a__A) = 1.0;                                         \
                                                                        \
} while(0)

// A_symm = symmetric(A)
// Matrix_symmetric(A,A_symm)
#define Matrix_symmetric(A,A_symm) do {                                 \
  symmetric_part((A_symm).m_pdata,(A).m_pdata,(A).m_row);               \
} while(0)
  
/* trA = tr_(A), trA = A_ii */
#define Matrix_trace(A,trA) do {                                        \
  long a__A;                                                            \
  trA = 0.0;                                                            \
  if(A.m_row != A.m_col || A.m_row ==0 || A.m_col ==0)                  \
    break;                                                              \
  for(a__A = 1; a__A <= A.m_row; a__A++)                                \
    trA += Mat_v(A, a__A, a__A);                                        \
} while(0)

#define Matrix_det(A, ddet) do {                                        \
                                                                        \
  if((A).m_row != (A).m_col || (A).m_row ==0 || (A).m_col ==0)          \
    break;                                                              \
                                                                        \
  if((A).m_row==1){                                                     \
    ddet = Mat_v(A, 1, 1);                                              \
    break;                                                              \
  }                                                                     \
                                                                        \
  if((A).m_row==2){                                                     \
    ddet = Mat_v(A, 1, 1)*Mat_v(A, 2, 2);                               \
    ddet -= Mat_v(A, 1, 2)*Mat_v(A, 2, 1);                              \
    break;                                                              \
  };                                                                    \
                                                                        \
  if((A).m_row==3){                                                     \
    ddet  = Mat_v(A, 1, 1)*Mat_v(A, 2, 2)*Mat_v(A, 3, 3);               \
    ddet += Mat_v(A, 1, 2)*Mat_v(A, 2, 3)*Mat_v(A, 3, 1);               \
    ddet += Mat_v(A, 1, 3)*Mat_v(A, 3, 2)*Mat_v(A, 2, 1);               \
    ddet -= Mat_v(A, 1, 3)*Mat_v(A, 2, 2)*Mat_v(A, 3, 1);               \
    ddet -= Mat_v(A, 1, 2)*Mat_v(A, 2, 1)*Mat_v(A, 3, 3);               \
    ddet -= Mat_v(A, 1, 1)*Mat_v(A, 3, 2)*Mat_v(A, 2, 3);               \
    break;                                                              \
  };                                                                    \
  if((A).m_row > 3){                                                    \
    printf("Matrix greater than [3x3] is not currently supported\n");   \
    break;                                                              \
  }                                                                     \
} while(0)                                                              

#define Matrix_inv_check_err(A, invA, info) do {                        \
  info = 0;                                                             \
  info += inverse((A).m_pdata,(A).m_row,(invA).m_pdata);                \
} while(0)

#define Matrix_inv(A, invA) do {                                        \
  inverse((A).m_pdata,(A).m_row,(invA).m_pdata);                        \
} while(0)

#define Matrix_inv_no_use(A, invA) do {                                 \
  double detA = 0.0;                                                    \
                                                                        \
  if((A).m_row != (A).m_col || (A).m_row ==0 || (A).m_col ==0)          \
    break;                                                              \
                                                                        \
  Matrix_check_null_and_redim(invA, (A).m_row, (A).m_col);              \
  if(A.m_row==1){                                                       \
    Mat_v(invA, 1, 1) = 1.0/Mat_v(A, 1, 1);                             \
    break;                                                              \
  }                                                                     \
                                                                        \
  Matrix_det(A, detA);                                                  \
  if(fabs(detA)<1.0e-15){                                               \
    printf("det(Matrix) = %f is close to zero\n", fabs(detA));          \
    break;                                                              \
  }                                                                     \
                                                                        \
  if((A).m_row==2){                                                     \
    Mat_v(invA, 1, 1) =  Mat_v(A, 2, 2)/detA;                           \
    Mat_v(invA, 2, 2) =  Mat_v(A, 1, 1)/detA;                           \
    Mat_v(invA, 1, 2) = -Mat_v(A, 1, 2)/detA;                           \
    Mat_v(invA, 2, 1) = -Mat_v(A, 2, 1)/detA;                           \
    break;                                                              \
  };                                                                    \
                                                                        \
  if((A).m_row==3){                                                     \
                                                                        \
    Mat_v(invA, 1, 1) = (Mat_v(A, 2, 2)*Mat_v(A, 3, 3)                  \
                       - Mat_v(A, 2, 3)*Mat_v(A, 3, 2))/detA;           \
    Mat_v(invA, 1, 2) = (Mat_v(A, 1, 3)*Mat_v(A, 3, 2)                  \
                       - Mat_v(A, 1, 2)*Mat_v(A, 3, 3))/detA;           \
    Mat_v(invA, 1, 3) = (Mat_v(A, 1, 2)*Mat_v(A, 2, 3)                  \
                       - Mat_v(A, 1, 3)*Mat_v(A, 2, 2))/detA;           \
    Mat_v(invA, 2, 1) = (Mat_v(A, 2, 3)*Mat_v(A, 3, 1)                  \
                       - Mat_v(A, 2, 1)*Mat_v(A, 3, 3))/detA;           \
    Mat_v(invA, 2, 2) = (Mat_v(A, 1, 1)*Mat_v(A, 3, 3)                  \
                       - Mat_v(A, 1, 3)*Mat_v(A, 3, 1))/detA;           \
    Mat_v(invA, 2, 3) = (Mat_v(A, 1, 3)*Mat_v(A, 2, 1)                  \
                       - Mat_v(A, 1, 1)*Mat_v(A, 2, 3))/detA;           \
    Mat_v(invA, 3, 1) = (Mat_v(A, 2, 1)*Mat_v(A, 3, 2)                  \
                       - Mat_v(A, 2, 2)*Mat_v(A, 3, 1))/detA;           \
    Mat_v(invA, 3, 2) = (Mat_v(A, 1, 2)*Mat_v(A, 3, 1)                  \
                       - Mat_v(A, 1, 1)*Mat_v(A, 3, 2))/detA;           \
    Mat_v(invA, 3, 3) = (Mat_v(A, 1, 1)*Mat_v(A, 2, 2)                  \
                       - Mat_v(A, 1, 2)*Mat_v(A, 2, 1))/detA;           \
    break;                                                              \
  };                                                                    \
  if((A).m_row > 3){                                                    \
    printf("Matrix greater than [3x3] is not currently supported\n");   \
    break;                                                              \
  }                                                                     \
} while(0)                                                              

#define Matrix_cleanup(A) do {                                          \
  (A).m_row = 0;                                                        \
  (A).m_col = 0;                                                        \
  if((A).m_pdata)                                                       \
    free((A).m_pdata);                                                  \
  (A).m_pdata = NULL;                                                   \
  if((A).temp)                                                          \
    free((A).temp);                                                     \
  (A).temp = NULL;                                                      \
} while(0)

#define Matrix_print(A) do {                                            \
  long a__A, b__B;                                                      \
  printf("[%ldx%ld] = \n", (A).m_row, (A).m_col);                       \
  for(a__A = 1; a__A <= (A).m_row; a__A++){                             \
    for(b__B = 1; b__B <= (A).m_col; b__B++){                           \
      printf("%e ", (double)Mat_v(A, a__A, b__B));                      \
    }                                                                   \
    printf("\n");                                                       \
  }                                                                     \
} while(0)

#define Matrix_print_name(A, name) do {                                 \
  long a__A, b__B;                                                      \
  printf("%s = [\n", name);                                             \
  for(a__A = 1; a__A <= (A).m_row; a__A++){                             \
    for(b__B = 1; b__B <= (A).m_col; b__B++){                           \
      printf("%e ", (double)Mat_v(A, a__A, b__B));                      \
    }                                                                   \
    if(a__A==(A).m_row)                                                 \
      printf("];\n");                                                   \
    else                                                                \
      printf("\n");                                                     \
  }                                                                     \
} while(0)

/* A <- B */
#define Matrix_copy(A,B) do {                                           \
    assert( ((A).m_row == (B).m_row) && ((A).m_col == (B).m_col) );     \
    if ( sizeof(*((A).m_pdata)) == sizeof(*((B).m_pdata)) ) {           \
      memcpy((A).m_pdata,(B).m_pdata,                                   \
             (A).m_row * (A).m_col * sizeof(*((A).m_pdata)));           \
    } else {                                                            \
      for ( size_t idx = 0, e = (A).m_row * (A).m_col;                  \
            idx < e; idx++) (A).m_pdata[idx] = (B).m_pdata[idx];        \
    }                                                                   \
} while (0)

/* A = bB */
#define Matrix_AeqB(A, b, B) do {                                       \
  long a__A;                                                            \
  long m_row = (B).m_row;                                               \
  long m_col = (B).m_col;                                               \
  Matrix_check_null_and_redim(A, (B).m_row, (B).m_col);                 \
                                                                        \
  for(a__A = 0; a__A < m_row*m_col; a__A++)                             \
    (A).m_pdata[a__A] = (B).m_pdata[a__A]*(b);                          \
} while(0)
  
/* A = transpose(A) */
#define Matrix_trans(A) do {                                            \
  long a__A, b__B;                                                      \
  long m_row = (A).m_row;                                               \
  long m_col = (A).m_col;                                               \
  (A).temp =  malloc((A).sizeof_T*m_row*m_col);                         \
                                                                        \
  for(a__A = 1; a__A <= m_row; a__A++){                                 \
    for(b__B = 1; b__B <= m_col; b__B++){                               \
      (A).temp[(b__B-1)*(A).m_col+(a__A-1)] = Mat_v(A, a__A, b__B);     \
    }                                                                   \
  }                                                                     \
  (A).m_row = m_col;                                                    \
  (A).m_col = m_row;                                                    \
  if((A).m_pdata)                                                       \
    free((A).m_pdata);                                                  \
  (A).m_pdata = (A).temp;                                               \
  (A).temp = NULL;                                                      \
} while(0)

/* A = bBT */
#define Matrix_AeqBT(A,b,B) do {                                        \
  long m_row = (B).m_row;                                               \
  long m_col = (B).m_col;                                               \
  Matrix_check_null_and_redim(A, m_col, m_row);                         \
                                                                        \
  for(int a__A = 1; a__A <= m_row; a__A++){                             \
    for(int b__B = 1; b__B <= m_col; b__B++){                           \
      Mat_v(A, b__B, a__A) = Mat_v(B,a__A, b__B)*(b);                   \
    }                                                                   \
  }                                                                     \
} while(0)

/* C = aA + bB */
#define Matrix_AplusB(C, a, A, b, B) do {                               \
  long i__I;                                                            \
  long m_row = (A).m_row;                                               \
  long m_col = (A).m_col;                                               \
  if((A).m_row != (B).m_row)                                            \
    break;                                                              \
  if((A).m_col != (B).m_col)                                            \
    break;                                                              \
                                                                        \
  Matrix_check_null_and_redim(C, (A).m_row, (A).m_col);                 \
  for(i__I = 0; i__I < m_row*m_col; i__I++)                             \
    (C).m_pdata[i__I] = (A).m_pdata[i__I]*(a) + (B).m_pdata[i__I]*(b);  \
} while(0)

#define Matrix_AOxB(C, A, B) do {                                       \
  for(int i__I=1; i__I<=3; i__I++)                                      \
  {                                                                     \
    for(int j__J=1; j__J<=3; j__J++)                                    \
    {                                                                   \
      for(int k__K=1; k__K<=3; k__K++)                                  \
      {                                                                 \
        for(int l__L=1; l__L<=3; l__L++)                                \
         Tns4_v(C, i__I,j__J,k__K,l__L) = Mat_v(A,i__I,j__J)*Mat_v(B,k__K,l__L);\
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while(0)

// C = A:B
#define Matrix_Tns4_dd_Tns2(C, A, B) do {                               \
  Matrix_init(C,0.0);                                                   \
  for(int k__K=1; k__K<=3; k__K++)                                      \
  {                                                                     \
    for(int l__L=1; l__L<=3; l__L++)                                    \
    {                                                                   \
      if(fabs(Mat_v(B,k__K,l__L))<1.0e-18)                              \
        continue;                                                       \
      for(int i__I=1; i__I<=3; i__I++)                                  \
      {                                                                 \
        for(int j__J=1; j__J<=3; j__J++)                                \
         Mat_v(C,i__I,j__J) += Tns4_v(A,i__I,j__J,k__K,l__L)*Mat_v(B,k__K,l__L);\
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while(0)

// C = A:B
#define Matrix_Tns2_dd_Tns4(C, A, B) do {                               \
  Matrix_init(C,0.0);                                                   \
  for(int i__I=1; i__I<=3; i__I++)                                      \
  {                                                                     \
    for(int j__J=1; j__J<=3; j__J++)                                    \
    {                                                                   \
      if(Mat_v(A,i__I,j__J)==0.0)                                       \
        continue;                                                       \
      for(int k__K=1; k__K<=3; k__K++)                                  \
      {                                                                 \
        for(int l__L=1; l__L<=3; l__L++)                                \
         Mat_v(C,k__K,l__L) += Mat_v(A,i__I,j__J)*Tns4_v(B,i__I,j__J,k__K,l__L);\
      }                                                                 \
    }                                                                   \
  }                                                                     \
} while(0)

#define Matrix_ddot(A,B,A_dd_B) do {                                    \
  A_dd_B = cblas_ddot((A).m_row*(A).m_col,(A).m_pdata,1,(B).m_pdata,1); \
} while(0)

#define Matrix_ddot_no_use(A,B,A_dd_B) do {                             \
  A_dd_B = 0.0;                                                         \
  for(int i__I = 0; i__I<(A).m_row*(A).m_col; i__I++)                   \
    A_dd_B += (A).m_pdata[i__I]*(B).m_pdata[i__I];                      \
} while(0)

// Ax=B, x = inv(A)*B
#define Matrix_solver(A,x,B) do {                                       \
  int *ipiv = malloc(sizeof(int)*(A).m_row);                            \
  int info;                                                             \
  Matrix_AeqB(x,1.0,B);                                                 \
  int N_N__ = (A).m_row;                                                \
  int n_n__ = (x).m_col;                                                \
  int m_m__ = (x).m_row;                                                \
  dgesv(&N_N__,&n_n__,(A).m_pdata,&N_N__,ipiv,(x).m_pdata,&m_m__,&info);\
  if(info >0)                                                           \
  {                                                                     \
    printf( "A is singular\n");                                         \
    printf( "the solution could not be computed.\n" );                  \
    exit( 1 );                                                          \
  }                                                                     \
  free(ipiv);                                                           \
} while(0)
          			      
/*C = aAxB + bC*/
// C[m,n] = a*A[m,k] x B[k,n] + b*C[m,n]
// cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,a,A,k,B,n,b,C,n);
#define Matrix_AxB(C, a, b, A, a__AT, B, b__BT) do {                                                                                                               \
  if(a__AT==0 && b__BT==0){                                                                                                                                        \
    Matrix_check_null_and_redim(C, (A).m_row, (B).m_col);                                                                                                    \
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(A).m_row,(B).m_col,(A).m_col,a,(A).m_pdata,(A).m_col,(B).m_pdata,(B).m_col,b,(C).m_pdata,(C).m_col);\
  }                                                                                                                                                          \
  if(a__AT==1 && b__BT==0){                                                                                                                                        \
    Matrix_check_null_and_redim(C, (A).m_col, (B).m_col);                                                                                                    \
    cblas_dgemm(CblasRowMajor,CblasTrans,  CblasNoTrans,(C).m_row,(C).m_col,(B).m_row,a,(A).m_pdata,(A).m_col,(B).m_pdata,(B).m_col,b,(C).m_pdata,(C).m_col);\
  }                                                                                                                                                          \
  if(a__AT==0 && b__BT==1){                                                                                                                                        \
    Matrix_check_null_and_redim(C, (A).m_row, (B).m_row);                                                                                                    \
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,  (C).m_row,(C).m_col,(A).m_col,a,(A).m_pdata,(A).m_col,(B).m_pdata,(B).m_col,b,(C).m_pdata,(C).m_col);\
  }                                                                                                                                                          \
  if(a__AT==1 && b__BT==1){                                                                                                                                        \
    Matrix_check_null_and_redim(C, (A).m_col, (B).m_row);                                                                                                    \
    cblas_dgemm(CblasRowMajor,CblasTrans,  CblasTrans,  (C).m_row,(C).m_col,(A).m_row,a,(A).m_pdata,(A).m_col,(B).m_pdata,(B).m_col,b,(C).m_pdata,(C).m_col);\
  }                                                                                                                                                          \
} while(0)

// D = a*A*B*C + b*D, a and b are scalar
#define Matrix_Tns2_AxBxC(D,a,b,A,B,C) do {                             \
                                                                        \
  Matrix_check_null_and_redim(D,3,3);                                   \
  Mat_v(D,1,1) = Mat_v(D,1,1)*(b) + (a)*(                                                                          \
               + Mat_v(C,1,1)*(Mat_v(A,1,1)*Mat_v(B,1,1) + Mat_v(A,1,2)*Mat_v(B,2,1) + Mat_v(A,1,3)*Mat_v(B,3,1))  \
               + Mat_v(C,2,1)*(Mat_v(A,1,1)*Mat_v(B,1,2) + Mat_v(A,1,2)*Mat_v(B,2,2) + Mat_v(A,1,3)*Mat_v(B,3,2))  \
               + Mat_v(C,3,1)*(Mat_v(A,1,1)*Mat_v(B,1,3) + Mat_v(A,1,2)*Mat_v(B,2,3) + Mat_v(A,1,3)*Mat_v(B,3,3)));\
  Mat_v(D,1,2) = Mat_v(D,1,2)*(b) + (a)*(                                                                          \
               + Mat_v(C,1,2)*(Mat_v(A,1,1)*Mat_v(B,1,1) + Mat_v(A,1,2)*Mat_v(B,2,1) + Mat_v(A,1,3)*Mat_v(B,3,1))  \
               + Mat_v(C,2,2)*(Mat_v(A,1,1)*Mat_v(B,1,2) + Mat_v(A,1,2)*Mat_v(B,2,2) + Mat_v(A,1,3)*Mat_v(B,3,2))  \
               + Mat_v(C,3,2)*(Mat_v(A,1,1)*Mat_v(B,1,3) + Mat_v(A,1,2)*Mat_v(B,2,3) + Mat_v(A,1,3)*Mat_v(B,3,3)));\
  Mat_v(D,1,3) = Mat_v(D,1,3)*(b) + (a)*(                                                                          \
               + Mat_v(C,1,3)*(Mat_v(A,1,1)*Mat_v(B,1,1) + Mat_v(A,1,2)*Mat_v(B,2,1) + Mat_v(A,1,3)*Mat_v(B,3,1))  \
               + Mat_v(C,2,3)*(Mat_v(A,1,1)*Mat_v(B,1,2) + Mat_v(A,1,2)*Mat_v(B,2,2) + Mat_v(A,1,3)*Mat_v(B,3,2))  \
               + Mat_v(C,3,3)*(Mat_v(A,1,1)*Mat_v(B,1,3) + Mat_v(A,1,2)*Mat_v(B,2,3) + Mat_v(A,1,3)*Mat_v(B,3,3)));\
  Mat_v(D,2,1) = Mat_v(D,2,1)*(b) + (a)*(                                                                          \
               + Mat_v(C,1,1)*(Mat_v(A,2,1)*Mat_v(B,1,1) + Mat_v(A,2,2)*Mat_v(B,2,1) + Mat_v(A,2,3)*Mat_v(B,3,1))  \
               + Mat_v(C,2,1)*(Mat_v(A,2,1)*Mat_v(B,1,2) + Mat_v(A,2,2)*Mat_v(B,2,2) + Mat_v(A,2,3)*Mat_v(B,3,2))  \
               + Mat_v(C,3,1)*(Mat_v(A,2,1)*Mat_v(B,1,3) + Mat_v(A,2,2)*Mat_v(B,2,3) + Mat_v(A,2,3)*Mat_v(B,3,3)));\
  Mat_v(D,2,2) = Mat_v(D,2,2)*(b) + (a)*(                                                                          \
               + Mat_v(C,1,2)*(Mat_v(A,2,1)*Mat_v(B,1,1) + Mat_v(A,2,2)*Mat_v(B,2,1) + Mat_v(A,2,3)*Mat_v(B,3,1))  \
               + Mat_v(C,2,2)*(Mat_v(A,2,1)*Mat_v(B,1,2) + Mat_v(A,2,2)*Mat_v(B,2,2) + Mat_v(A,2,3)*Mat_v(B,3,2))  \
               + Mat_v(C,3,2)*(Mat_v(A,2,1)*Mat_v(B,1,3) + Mat_v(A,2,2)*Mat_v(B,2,3) + Mat_v(A,2,3)*Mat_v(B,3,3)));\
  Mat_v(D,2,3) = Mat_v(D,2,3)*(b) + (a)*(                                                                          \
               + Mat_v(C,1,3)*(Mat_v(A,2,1)*Mat_v(B,1,1) + Mat_v(A,2,2)*Mat_v(B,2,1) + Mat_v(A,2,3)*Mat_v(B,3,1))  \
               + Mat_v(C,2,3)*(Mat_v(A,2,1)*Mat_v(B,1,2) + Mat_v(A,2,2)*Mat_v(B,2,2) + Mat_v(A,2,3)*Mat_v(B,3,2))  \
               + Mat_v(C,3,3)*(Mat_v(A,2,1)*Mat_v(B,1,3) + Mat_v(A,2,2)*Mat_v(B,2,3) + Mat_v(A,2,3)*Mat_v(B,3,3)));\
  Mat_v(D,3,1) = Mat_v(D,3,1)*(b) + (a)*(                                                                          \
               + Mat_v(C,1,1)*(Mat_v(A,3,1)*Mat_v(B,1,1) + Mat_v(A,3,2)*Mat_v(B,2,1) + Mat_v(A,3,3)*Mat_v(B,3,1))  \
               + Mat_v(C,2,1)*(Mat_v(A,3,1)*Mat_v(B,1,2) + Mat_v(A,3,2)*Mat_v(B,2,2) + Mat_v(A,3,3)*Mat_v(B,3,2))  \
               + Mat_v(C,3,1)*(Mat_v(A,3,1)*Mat_v(B,1,3) + Mat_v(A,3,2)*Mat_v(B,2,3) + Mat_v(A,3,3)*Mat_v(B,3,3)));\
  Mat_v(D,3,2) = Mat_v(D,3,2)*(b) + (a)*(                                                                          \
               + Mat_v(C,1,2)*(Mat_v(A,3,1)*Mat_v(B,1,1) + Mat_v(A,3,2)*Mat_v(B,2,1) + Mat_v(A,3,3)*Mat_v(B,3,1))  \
               + Mat_v(C,2,2)*(Mat_v(A,3,1)*Mat_v(B,1,2) + Mat_v(A,3,2)*Mat_v(B,2,2) + Mat_v(A,3,3)*Mat_v(B,3,2))  \
               + Mat_v(C,3,2)*(Mat_v(A,3,1)*Mat_v(B,1,3) + Mat_v(A,3,2)*Mat_v(B,2,3) + Mat_v(A,3,3)*Mat_v(B,3,3)));\
  Mat_v(D,3,3) = Mat_v(D,3,3)*(b) + (a)*(                                                                          \
               + Mat_v(C,1,3)*(Mat_v(A,3,1)*Mat_v(B,1,1) + Mat_v(A,3,2)*Mat_v(B,2,1) + Mat_v(A,3,3)*Mat_v(B,3,1))  \
               + Mat_v(C,2,3)*(Mat_v(A,3,1)*Mat_v(B,1,2) + Mat_v(A,3,2)*Mat_v(B,2,2) + Mat_v(A,3,3)*Mat_v(B,3,2))  \
               + Mat_v(C,3,3)*(Mat_v(A,3,1)*Mat_v(B,1,3) + Mat_v(A,3,2)*Mat_v(B,2,3) + Mat_v(A,3,3)*Mat_v(B,3,3)));\
} while(0)

#endif
