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

/*==== Tangent of Material stiffness functions ====*/
void SEDF_d3W_dC3_Mooney_Rivlin(double *C,
                                 MATERIAL_ELASTICITY const *mat,
                                 double *K)
{
  Matrix(double) temp, invC;
  Matrix_construct_redim(double, invC, DIM_3, DIM_3);

  temp.m_row = temp.m_col = DIM_3;
  temp.m_pdata = C;
  
  Matrix_inv(temp, invC);

  double detC, trC, CC;
  Matrix_det(temp, detC);
  Matrix_trace(temp,trC);
  Matrix_ddot(temp,temp,CC);  

  const double *C_I = invC.m_pdata; /* get constatnt pointer */
  const double J23 = pow(detC,-1.0/3.0);
  const double J43 = J23*J23;

  for(int i=0; i<DIM_3; i++)
  {
    for(int j=0; j<DIM_3; j++)
    {
      const int ij = DIM_3*i + j;
      for(int k=0; k<DIM_3; k++)
      {
        const int ik = DIM_3*i + k;
        for(int l=0; l<DIM_3; l++)
        {
          const int kl = DIM_3*k + l;
          const int lj = DIM_3*l + j;
          for(int r=0; r<DIM_3; r++)
          {
            const int ir = DIM_3*i + r;
            const int kr = DIM_3*k + r;
            const int lr = DIM_3*l + r;
            for(int s=0; s<DIM_3; s++)
            {
              const int rs = DIM_3*r + s;
              const int sj = DIM_3*s + j;
              const int sk = DIM_3*s + k;
              const int sl = DIM_3*s + l;              
              const int ijklrs = DIM_3x3x3x3*3*i + DIM_3x3x3x3*j + DIM_3x3x3*k + DIM_3x3*l + DIM_3*r + s;
              
              /*============ mu_10 term ==================*/
              double M_10_term = 0.0;
              M_10_term += 1./9.*((r==s)*C_I[ij] - trC*C[ir]*C_I[sj])*C_I[kl];
              
              M_10_term -= (1./9.*trC*C_I[ij]-1./3.*(i==j))*C_I[kr]*C_I[sl];
              
              M_10_term += 1./3.*((r==s)*C_I[ik]*C_I[lj]
                      - trC*C_I[ir]*C_I[sk]*C_I[lj]
                      - trC*C_I[ik]*C_I[lr]*C_I[sj]
                      + (k==l)*C_I[ir]*C_I[sj]);
              
              M_10_term -= (1./27.*trC*C_I[ij]-1./9.*(i==j))*C_I[kl]*C_I[rs];
              
              M_10_term -= 1./9.*(trC*C_I[ik]*C_I[lj]-(k==l)*C_I[ij])*C_I[rs];
              
              M_10_term *= 8.*mat->m10*J23;
              
              /*============ mu_01 term ==================*/
              double M_01_term = 0.0;
              M_01_term += 2./3.*((i==r)*(s==j) - (r==s)*(i==j))*C_I[kl];
              
              M_01_term -= 2./3.*(C[ij]-trC*(i==j))*C_I[kr]*C_I[sl];
              
              M_01_term -= 2.*(C[rs]-trC*(r==s))*(2./9.*C_I[ij]*C_I[kl]
                      + 1./3.*C_I[ik]*C_I[lj]);
              
              M_01_term += (CC-trC*trC)*(2./9.*C_I[ir]*C_I[sj]*C_I[kl]
                      + 2./9.*C_I[ij]*C_I[kr]*C_I[sl]
                      + 1./3.*C_I[ir]*C_I[sk]*C_I[lj]
                      + 1./3.*C_I[ik]*C_I[lr]*C_I[sj]);
              
              M_01_term += 2./3.*((k==r)*(s==l) - (r==s)*(k==l))*C_I[ij];
              
              M_01_term -= 2./3.*(C[kl]-trC*(k==l))*C_I[ir]*C_I[sj];
              
              M_01_term -= 2./3.*((i==j)*(k==l)-(i==k)*(l==j))*C_I[rs];
              
              M_01_term -= 4./9.*(C[ij]-trC*(i==j))*C_I[kl]*C_I[rs];
              
              M_01_term += 2./3.*(CC-trC*trC)*(2./9.*C_I[ij]*C_I[kl]
                      + 1./3.*C_I[ik]*C_I[lj])*C_I[rs];
              
              M_01_term -= 4./9.*(C[kl]-trC*(k==l))*C_I[ij]*C_I[rs];
              
              M_01_term *= 8.*mat->m01*J43;
              
              /*============ K ==================*/
              K[ijklrs] = M_10_term + M_01_term;
            }
          }
        }
      }
    }
  }
  Matrix_cleanup(invC);
}

void SEDF_d3W_dC3_Linear(double *C,
		                     MATERIAL_ELASTICITY const *mat,
		                     double *K)
{
  memset(K,0,729*sizeof(double));
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

/*==== d3UdJ3 functions ===*/
void SEDF_d3UdJ3_Common_new(double *d3UdJ3, double J)
{
  (*d3UdJ3) = 0.0;
}

void SEDF_d3UdJ3_Doll_Schweizerhof_7(double *d3UdJ3, double J)
{
  (*d3UdJ3) = (exp(-1. + J) - 2.0/(J*J*J))*0.5;
}

void SEDF_d3UdJ3_Doll_Schweizerhof_8(double *d3UdJ3, double J)
{
  (*d3UdJ3) = -0.5*(J + 2.0)/(J*J*J) ;
}
