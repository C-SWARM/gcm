#include "constitutive_model.h"
#include "strain_energy_density_function.h"
#include "material_properties.h"

/*==== Potential energy ====*/
void SEDF_devPotential_Mooney_Rivlin(double *C_in,
			     MATERIAL_ELASTICITY const *mat,
			     double *W)
{
  Matrix(double) C, Ch;
  Matrix_construct_redim(double, Ch, DIM_3, DIM_3);

  C.m_row = C.m_col = DIM_3;
  C.m_pdata = C_in;
  
  double detC, trCh, ChCh;
  Matrix_det(C, detC);
  
  Matrix_AeqB(Ch,pow(detC,-1.0/3.0), C);
  Matrix_trace(Ch,trCh);
  Matrix_ddot(Ch,Ch,ChCh);
  
  *W = mat->m10*(trCh-3.0) + 0.5*mat->m01*(trCh*trCh - ChCh-6.0);
  Matrix_cleanup(Ch);
}	


void SEDF_devPotential_Linear(double *C_in,
			 MATERIAL_ELASTICITY const *mat,
			 double *W)
{
  /* W_hat = int(S_hat)dC || W_hat = 0|C=I */
  /* S_hat = 2G E || E = 1/2 (C-1) */

  Matrix(double) C, Ch;

  C.m_row = C.m_col = DIM_3;
  C.m_pdata = C_in;
  
  double CC, trC;
  Matrix_ddot(C,C,CC);
  Matrix_trace(C,trC);

   *W = 0.25*mat->G*(CC - 2.0*trC + 3.0);
}		     
			     
/*==== Deviatoric stress functions ====*/
void SEDF_devStress_Mooney_Rivlin(double *C_in,
			     MATERIAL_ELASTICITY const *mat,
			     double *S)
{
  Matrix(double) C, ident, invC;
  Matrix_construct_redim(double, ident, DIM_3, DIM_3);
  Matrix_construct_redim(double, invC, DIM_3, DIM_3);

  Matrix_eye(ident,DIM_3);

  C.m_row = C.m_col = DIM_3;
  C.m_pdata = C_in;
  
  Matrix_inv(C, invC);
  double detC, trC, CC;
  Matrix_det(C, detC);
  Matrix_trace(C,trC);
  Matrix_ddot(C,C,CC);

  for(int i=0; i<DIM_3x3; i++){
    S[i] =( 
	   (2.0*mat->m10*pow(detC,-1./3.)
	    *(ident.m_pdata[i] - 1./3.*trC*invC.m_pdata[i]))	    
	   +(2.0*mat->m01*pow(detC,-2./3.)
	     *(trC*ident.m_pdata[i] - C.m_pdata[i]
	       + 1./3.*(CC-trC*trC)*invC.m_pdata[i]))
	    );
  }

 Matrix_cleanup(invC);
 Matrix_cleanup(ident);
}

void SEDF_devStress_Linear(double *C,
		      MATERIAL_ELASTICITY const *mat,
		      double *S)
{
  double ident[DIM_3x3];
  ident[1] = ident[2] = ident[3] = 0.0;
  ident[5] = ident[6] = ident[7] = 0.0;  
  ident[0] = ident[4] = ident[8] = 1.0;

  for(int i=0; i<9; i++){
    S[i] = mat->mu*(C[i]-ident[i]);
  }
}

/*==== Material stiffness functions ====*/
void SEDF_matStiff_Mooney_Rivlin(double *C_in,
			    MATERIAL_ELASTICITY const *mat,
			    double *L_out)
{
  Matrix(double) C, L, ident, invC;
  Matrix_construct_redim(double, ident, DIM_3, DIM_3);
  Matrix_construct_redim(double, invC, DIM_3, DIM_3);

  Matrix_eye(ident,DIM_3);

  C.m_row = C.m_col = DIM_3;
  C.m_pdata = C_in;
  
  L.m_row = DIM_3x3x3x3; L.m_col = 1;
  L.m_pdata = L_out;  
  
  Matrix_inv(C, invC);

  double detC, trC, CC;
  Matrix_det(C, detC);
  Matrix_trace(C,trC);
  Matrix_ddot(C,C,CC);  

  double detC1,detC2;
  detC1 = pow(detC,-1./3.);
  detC2 = detC1*detC1;

  for(int i=1; i<=DIM_3; i++)
  {
    for(int j=1; j<=DIM_3; j++)
    {
      for(int k=1; k<=DIM_3; k++)
      {
	      for(int l=1; l<=DIM_3; l++)
	      {
	        Tns4_v(L,i,j,k,l) = (4.0*mat->m10*detC1
			       * (
			      	  (1./9.*trC*Mat_v(invC,i,j)- 1./3.*Mat_v(ident,i,j))*Mat_v(invC,k,l)				 
			      	  + 1./3.*(trC*Mat_v(invC,i,k)*Mat_v(invC,l,j)-Mat_v(ident,k,l)*Mat_v(invC,i,j))
			      	 )
			       +(4.0*mat->m01*detC2*(
			      	    Mat_v(ident,i,j)*Mat_v(ident,k,l)-Mat_v(ident,i,k)*Mat_v(ident,l,j)				    
			      	    + 2./3.*((Mat_v(C,i,j)-trC*Mat_v(ident,i,j))*Mat_v(invC,k,l))
			      	    - ((CC-trC*trC)*(2./9.*Mat_v(invC,i,j)*Mat_v(invC,k,l)+1./3.*Mat_v(invC,i,k)*Mat_v(invC,l,j)))
			      	    + 2./3.*((Mat_v(C,k,l)-trC*Mat_v(ident,k,l))*Mat_v(invC,i,j)))
			      	 ));
        }
      }
    }
  }
  Matrix_cleanup(invC);
  Matrix_cleanup(ident);
}

void SEDF_matStiff_Linear(double *C,
		     MATERIAL_ELASTICITY const *mat,
		     double *L)
{
  double ident[DIM_3x3];
  ident[1] = ident[2] = ident[3] = 0.0;
  ident[5] = ident[6] = ident[7] = 0.0;  
  ident[0] = ident[4] = ident[8] = 1.0;
  
  for(int i=0; i<DIM_3; i++){
    for(int j=0; j<DIM_3; j++){
      for(int k=0; k<DIM_3; k++){
	      for(int l=0; l<DIM_3; l++){
	  int id4 = i*DIM_3x3x3+j*DIM_3x3+k*DIM_3+l;
	  L[id4] = (2*mat->mu*ident[i*DIM_3+k]
			       *ident[j*DIM_3+l]);
	}}}}
}

/*==== Volumetric Potetntial Functions ====*/
void SEDF_U_Common(double *U, double const J)
{
  *U = (0.5*(J-1.0)*(J-1.0));
}

void SEDF_U_Doll_Schweizerhof_7(double *U, double const J)
{
  *U = (0.5*(exp(J-1.0)-log(J)-1.0));
}

void SEDF_U_Doll_Schweizerhof_8(double *U, double const J)
{
  *U = (0.5*(J-1.0)*log(J));
}

/*==== dUdJ functions ====*/
void SEDF_dUdJ_Common(double *dUdJ, double J)
{
  (*dUdJ) = J-1.0;
}

void SEDF_dUdJ_Doll_Schweizerhof_7(double *dUdJ, double J)
{
  (*dUdJ) = (exp(-1.0 + J) - 1.0/J)*0.5;
}

void SEDF_dUdJ_Doll_Schweizerhof_8(double *dUdJ, double J)
{
  (*dUdJ) = (-1.0 + J)/(2.0*J) + log(J)*0.5;
}

/*==== d2UdJ2 functions ===*/
void SEDF_d2UdJ2_Common_new(double *d2UdJ2, double J)
{
  (*d2UdJ2) = 1.0;
}

void SEDF_d2UdJ2_Doll_Schweizerhof_7(double *d2UdJ2, double J)
{
  (*d2UdJ2) = (exp(-1.0 + J) + 1./(J*J))*0.5;
}

void SEDF_d2UdJ2_Doll_Schweizerhof_8(double *d2UdJ2, double J)
{
  (*d2UdJ2) = 1.0/J - (-1.0 + J)/(2.0*J*J);
}


