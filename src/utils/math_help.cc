#include "math_help.h"
#include "mkl_lapack.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
 
double det2x2(const double *mat)
{
  return (mat[0]*mat[3] - mat[1]*mat[2]);
}

double det3x3(const double *mat)
{
  return ( ( mat[0]*(mat[4]*mat[8] - mat[5]*mat[7]) )
	   + mat[1]*(mat[5]*mat[6] - mat[3]*mat[8])
	   + mat[2]*(mat[3]*mat[7] - mat[4]*mat[6]) );
}

double det4x4(const double *mat)
{
  return (mat[1]*mat[11]*mat[14]*mat[4]
	  - mat[1]*mat[10]*mat[15]*mat[4]
	  - mat[11]*mat[13]*mat[2]*mat[4]
	  + mat[10]*mat[13]*mat[3]*mat[4]
	  - mat[0]*mat[11]*mat[14]*mat[5]
	  + mat[0]*mat[10]*mat[15]*mat[5]
	  + mat[11]*mat[12]*mat[2]*mat[5]
	  - mat[10]*mat[12]*mat[3]*mat[5]
	  - mat[1]*mat[11]*mat[12]*mat[6]
	  + mat[0]*mat[11]*mat[13]*mat[6]
	  + mat[1]*mat[10]*mat[12]*mat[7]
	  - mat[0]*mat[10]*mat[13]*mat[7]
	  - mat[15]*mat[2]*mat[5]*mat[8]
	  + mat[14]*mat[3]*mat[5]*mat[8]
	  + mat[1]*mat[15]*mat[6]*mat[8]
	  - mat[13]*mat[3]*mat[6]*mat[8]
	  - mat[1]*mat[14]*mat[7]*mat[8]
	  + mat[13]*mat[2]*mat[7]*mat[8]
	  + mat[15]*mat[2]*mat[4]*mat[9]
	  - mat[14]*mat[3]*mat[4]*mat[9]
	  - mat[0]*mat[15]*mat[6]*mat[9]
	  + mat[12]*mat[3]*mat[6]*mat[9]
	  + mat[0]*mat[14]*mat[7]*mat[9]
	  - mat[12]*mat[2]*mat[7]*mat[9]);
}


int inv2x2(const double *mat,
	   double *mat_inv)
{
  double A = det2x2(mat);
  if(A != 0.0){
    mat_inv[0] = 1.0/A*mat[3];
    mat_inv[1] = -1.0/A*mat[1];
    mat_inv[2] = -1.0/A*mat[2];
    mat_inv[3] = 1.0/A*mat[0];
    return 0;
  } else return 1;
}

int inv3x3(const double *mat,
	   double *mat_inv)
{
  double A = det3x3(mat);

  if (fabs(A) >= 1.0e-15){
    /* source: Wikipedia */
    mat_inv[0] = 1/A*(mat[4]*mat[8] - mat[5]*mat[7]);
    mat_inv[1] = -1/A*(mat[1]*mat[8] - mat[2]*mat[7]);
    mat_inv[2] = 1/A*(mat[1]*mat[5] - mat[2]*mat[4]);

    mat_inv[3] = -1/A*(mat[3]*mat[8] - mat[5]*mat[6]);
    mat_inv[4] = 1/A*(mat[0]*mat[8] - mat[2]*mat[6]);
    mat_inv[5] = -1/A*(mat[0]*mat[5] - mat[2]*mat[3]);

    mat_inv[6] = 1/A*(mat[3]*mat[7] - mat[4]*mat[6]);
    mat_inv[7] = -1/A*(mat[0]*mat[7] - mat[1]*mat[6]);
    mat_inv[8] = 1/A*(mat[0]*mat[4] - mat[1]*mat[3]);
    return 0;
  } else return 1;
}

int inv4x4(const double *mat,
	   double *mat_inv)
{
  double A = det4x4(mat);

  if(A != 0.0){
  mat_inv[0] = (-mat[11]*mat[14]*mat[5]
		+ mat[10]*mat[15]*mat[5]
		+ mat[11]*mat[13]*mat[6]
		- mat[10]*mat[13]*mat[7]
		- mat[15]*mat[6]*mat[9]
		+ mat[14]*mat[7]*mat[9]);

  mat_inv[1] = (mat[1]*mat[11]*mat[14]
		- mat[1]*mat[10]*mat[15]
		- mat[11]*mat[13]*mat[2]
		+ mat[10]*mat[13]*mat[3]
		+ mat[15]*mat[2]*mat[9]
		- mat[14]*mat[3]*mat[9]);

  mat_inv[2] = (-mat[15]*mat[2]*mat[5]
		+ mat[14]*mat[3]*mat[5]
		+ mat[1]*mat[15]*mat[6]
		- mat[13]*mat[3]*mat[6]
		- mat[1]*mat[14]*mat[7]
		+ mat[13]*mat[2]*mat[7]);

  mat_inv[3] = (mat[11]*mat[2]*mat[5]
		- mat[10]*mat[3]*mat[5]
		- mat[1]*mat[11]*mat[6]
		+ mat[1]*mat[10]*mat[7]
		+ mat[3]*mat[6]*mat[9]
		- mat[2]*mat[7]*mat[9]);

  mat_inv[4] = (mat[11]*mat[14]*mat[4]
		- mat[10]*mat[15]*mat[4]
		- mat[11]*mat[12]*mat[6]
		+ mat[10]*mat[12]*mat[7]
		+ mat[15]*mat[6]*mat[8]
		- mat[14]*mat[7]*mat[8]);

  mat_inv[5] = (-mat[0]*mat[11]*mat[14]
		+ mat[0]*mat[10]*mat[15]
		+ mat[11]*mat[12]*mat[2]
		- mat[10]*mat[12]*mat[3]
		- mat[15]*mat[2]*mat[8]
		+ mat[14]*mat[3]*mat[8]);

  mat_inv[6] = (mat[15]*mat[2]*mat[4]
		- mat[14]*mat[3]*mat[4]
		- mat[0]*mat[15]*mat[6]
		+ mat[12]*mat[3]*mat[6]
		+ mat[0]*mat[14]*mat[7]
		- mat[12]*mat[2]*mat[7]);

  mat_inv[7] = (-mat[11]*mat[2]*mat[4]
		+ mat[10]*mat[3]*mat[4]
		+ mat[0]*mat[11]*mat[6]
		- mat[0]*mat[10]*mat[7]
		- mat[3]*mat[6]*mat[8]
		+ mat[2]*mat[7]*mat[8]);

  mat_inv[8] = (-mat[11]*mat[13]*mat[4]
		+ mat[11]*mat[12]*mat[5]
		- mat[15]*mat[5]*mat[8]
		+ mat[13]*mat[7]*mat[8]
		+ mat[15]*mat[4]*mat[9]
		- mat[12]*mat[7]*mat[9]);

  mat_inv[9] = (-mat[1]*mat[11]*mat[12]
		+ mat[0]*mat[11]*mat[13]
		+ mat[1]*mat[15]*mat[8]
		- mat[13]*mat[3]*mat[8]
		- mat[0]*mat[15]*mat[9]
		+ mat[12]*mat[3]*mat[9]);

  mat_inv[10] = (-mat[1]*mat[15]*mat[4]
		 + mat[13]*mat[3]*mat[4]
		 + mat[0]*mat[15]*mat[5]
		 - mat[12]*mat[3]*mat[5]
		 + mat[1]*mat[12]*mat[7]
		 - mat[0]*mat[13]*mat[7]);

  mat_inv[11] = (mat[1]*mat[11]*mat[4]
		 - mat[0]*mat[11]*mat[5]
		 + mat[3]*mat[5]*mat[8]
		 - mat[1]*mat[7]*mat[8]
		 - mat[3]*mat[4]*mat[9]
		 + mat[0]*mat[7]*mat[9]);

  mat_inv[12] = (mat[10]*mat[13]*mat[4]
		 - mat[10]*mat[12]*mat[5]
		 + mat[14]*mat[5]*mat[8]
		 - mat[13]*mat[6]*mat[8]
		 - mat[14]*mat[4]*mat[9]
		 + mat[12]*mat[6]*mat[9]);

  mat_inv[13] = (mat[1]*mat[10]*mat[12]
		 - mat[0]*mat[10]*mat[13]
		 - mat[1]*mat[14]*mat[8]
		 + mat[13]*mat[2]*mat[8]
		 + mat[0]*mat[14]*mat[9]
		 - mat[12]*mat[2]*mat[9]);

  mat_inv[14] = (mat[1]*mat[14]*mat[4]
		 - mat[13]*mat[2]*mat[4]
		 - mat[0]*mat[14]*mat[5]
		 + mat[12]*mat[2]*mat[5]
		 - mat[1]*mat[12]*mat[6]
		 + mat[0]*mat[13]*mat[6]);

  mat_inv[15] = (-mat[1]*mat[10]*mat[4]
		 + mat[0]*mat[10]*mat[5]
		 - mat[2]*mat[5]*mat[8]
		 + mat[1]*mat[6]*mat[8]
		 + mat[2]*mat[4]*mat[9]
		 - mat[0]*mat[6]*mat[9]);
  for(int i=0; i<16; i++) mat_inv[i] /= A;
  return 0;
  } else return 1;
}

int inverse(double const* A,
	    const int M,
	    double *A_I)
{
  if (M <= 0){return 1;}
  int info;
  int *iPerm;
  int lwork = M*M;
  double *work;

  info = 0;
  switch(M){
  case 1:
    if(A[0] == 0.0){
      info = 1;
      break;
    } else{
      A_I[0] = 1.0/A[0];
      break;
    }

  case 2: info = inv2x2(A,A_I); break;
  case 3: info = inv3x3(A,A_I); break;
  case 4: info = inv4x4(A,A_I); break;
  default:
    iPerm = (int *) malloc(M*sizeof(int));
    work = (double *) malloc(lwork*sizeof(double));

    memcpy(A_I,A,M*M*sizeof(double));

    /* Factor into U matrix */
#ifdef ARCH_BGQ
    dgetrf(M,M,A_I,M,iPerm,&info);
#else
    dgetrf(&M,&M,A_I,&M,iPerm,&info);
#endif
    
    if(DEBUG_PRINT_ERROR)
    {
      if(info<0){
        printf("WARNING: illegal parameter given"
	        " to dgetrf at position %d.\n",info);
      } else if(info>0){
        printf("WARNING: factor U is singular.\n");
      }
    }

    /* Compute inverse using factored matrix */
#ifdef ARCH_BGQ
    dgetri(M,A_I,M,iPerm,work,lwork,&info);
#else
    dgetri(&M,A_I,&M,iPerm,work,&lwork,&info);
#endif

    if(DEBUG_PRINT_ERROR)
    {  
      if(info<0){
        printf("WARNING: illegal parameter given"
	        " to dgetri at position %d.\n",info);
	    }    
    } 
  
    free(iPerm);
    free(work);
    break;
  }/* switch M */

  if(DEBUG_PRINT_ERROR)
  {
    if(info != 0){
      if(info > 0){
        printf("ERROR: Matrix is singular, inverse not computed.\n");
      } else {
        printf("ERROR: Error (%d) in inverse routine.\n",info);
      }
    }
  }

  return info;
}

void symmetric_part(double *sym,
		    const double *mat,
		    const int dim)
{
  for(int i=0; i<dim; i++){
    for(int j=0; j<dim; j++){
      sym[i*dim+j] = 0.5*(mat[i*dim+j]
					 + mat[j*dim+i]);
    }
  }
}