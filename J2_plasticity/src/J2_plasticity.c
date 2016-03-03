#include "constitutive_model.h"
#include "J2_plasticity.h"
#include "hyperelasticity.h"

#define DIM_3        3
#define DIM_3x3      9
#define DIM_3x3x3   27
#define DIM_3x3x3x3 81
#define J2P_INT_ALG_TOL 1.0e-10

#define idx_2(row,col) ((unsigned int) ((row)*3+(col)))

int compute_push_forward(Matrix(double) *F,
                         Matrix(double) *A,
                         Matrix(double) *a)
{  
  int err = 0;
  Matrix(double) wkspc;
  Matrix_construct_redim(double,wkspc,DIM_3,DIM_3);
  Matrix_AxB(wkspc,1.0,0.0,*F,0,*A,0);
  Matrix_AxB(*a,1.0,0.0,wkspc,0,*F,1);    
  Matrix_cleanup(wkspc);
  return err;
}

int compute_pull_back(Matrix(double) *F,
                      Matrix(double) *a,
                      Matrix(double) *A)
{
  int err = 0;
  Matrix(double) FI;
  Matrix_construct_init(double,FI,DIM_3,DIM_3,0.0);
  Matrix_inv(*F,FI);
  compute_push_forward(&FI, a, A);
  Matrix_cleanup(FI);
  return err;
}

int compute_pull_Tensor4(Matrix(double) *FI,
                         Matrix(double) *aep,
                         Matrix(double) *Aep)
{
  int err = 0;
  
  Matrix_init(*Aep,0.0);

  for (int i = 1; i <= DIM_3; i++) {
    for (int j = 1; j <= DIM_3; j++) {
      for (int k = 1; k <= DIM_3; k++) {
        for (int l = 1; l <= DIM_3; l++) {
          for (int m = 1; m <= DIM_3; m++) {
            for (int n = 1; n <= DIM_3; n++) {
              for (int o = 1; o <= DIM_3; o++) {
                for (int p = 1; p <= DIM_3; p++) {
                  Tns4_v(*Aep,i,j,k,l) += (Mat_v(*FI, i,m)*Mat_v(*FI, j,n)
                                * Mat_v(*FI, k,o) * Mat_v(*FI, l,p)*Tns4_v(*aep,m,n,o,p));
                }
              }
            }
          }
        }
      }
    }
  }

  return err;
}

// dev(a) = a - 1/3 tr(a) i
int compute_dev(Matrix(double) *a,
                Matrix(double) *dev_a)
{
  int err = 0;
  double tra = (a->m_pdata[0] + a->m_pdata[4] + a->m_pdata[8])/3.0;
  Matrix_AeqB(*dev_a,1.0,*a);

  dev_a->m_pdata[0] -= tra;
  dev_a->m_pdata[4] -= tra;
  dev_a->m_pdata[8] -= tra;
  return err;
}

int compute_s0(double G,
               Matrix(double) *bbar,
               Matrix(double) *s0)
{
  int err = 0;
  double tra = (bbar->m_pdata[0] + bbar->m_pdata[4] + bbar->m_pdata[8])/3.0;
  Matrix_AeqB(*s0,G,*bbar);

  s0->m_pdata[0] -= G*tra;
  s0->m_pdata[4] -= G*tra;
  s0->m_pdata[8] -= G*tra;  
  return err;
}

/* bbar = J^-(2/3) F F' */
int compute_bbar_of_J(Matrix(double) *bbar, Matrix(double) *F, double J)
{
  int err = 0;
  double J23 = pow(J,-2.0/3.0);  
  Matrix_AxB(*bbar,J23,0.0,*F,0,*F,1);
  return err;
}

int compute_bbar(Matrix(double) *bbar, Matrix(double) *F)
{
  int err = 0;
  
  double J = 0.0;
  Matrix_det(*F, J);
  double J23 = pow(J,-2.0/3.0);
  
  Matrix_AxB(*bbar,J23,0.0,*F,0,*F,1);
  return err;
}

int compute_Fubar(Matrix(double) *F,
                  Matrix(double) *Fn,
                  Matrix(double) *Fubar)
{
  int err = 0;
  Matrix(double) FnI;
  Matrix_construct_redim(double,FnI,DIM_3,DIM_3);
  Matrix_inv(*Fn,FnI);
  Matrix_AxB(*Fubar,1.0,0.0,*F,0,FnI,0);
  double J = 0.0;
  Matrix_det(*Fubar,J);
  double Ju13 = pow(J, -1.0/3.0);
  for(int i = 0; i < DIM_3x3; i++) 
    Fubar->m_pdata[i] *= Ju13;

  Matrix_cleanup(FnI);
  return err;
}

int compute_sp_tr(Matrix(double) *F,
                  Matrix(double) *Fn,
                  Matrix(double) *spn,
                  Matrix(double) *sp_tr)
{
  int err = 0;
  Matrix(double) Fubar;
  Matrix_construct_redim(double,Fubar,DIM_3,DIM_3);
    
  err += compute_Fubar(F,Fn,&Fubar);
  err += compute_push_forward(&Fubar,spn,sp_tr);
  
  double tr = (sp_tr->m_pdata[0] + sp_tr->m_pdata[4] + sp_tr->m_pdata[8])/3.0;
  sp_tr->m_pdata[0] -= tr;
  sp_tr->m_pdata[4] -= tr;
  sp_tr->m_pdata[8] -= tr;

  Matrix_cleanup(Fubar);
  return err;
}

double compute_normal(MATERIAL_J2_PLASTICITY *J2P,
                      MATERIAL_ELASTICITY *mat_e,
                      Matrix(double) *s_tr, 
                      Matrix(double) *sp_tr, 
                      Matrix(double) *n)
{
  double G = mat_e->G;  
  double coef = (J2P->hp)/(3.0*G)*(1.0 - J2P->beta);
  double nrm = 0;

  /* compute ksi_tr and ||ksi_tr|| simultaneously */
  for(int i = 0; i < DIM_3x3; i++)
  {
    n->m_pdata[i] = s_tr->m_pdata[i] - coef*(sp_tr->m_pdata[i]);
    nrm += (n->m_pdata[i])*(n->m_pdata[i]);
  }
  nrm = sqrt(nrm);

  /* compute normal: n = ksi_tr / ||ksi_tr|| */
  if(nrm > 0) 
  {
    for(int i = 0; i < DIM_3x3; i++) n->m_pdata[i] /= nrm;
  }
  return nrm;  
};

double phi_yield_functioncompute_normal(MATERIAL_J2_PLASTICITY *J2P,
                                        MATERIAL_ELASTICITY *mat_e,
                                        Matrix(double) *s_tr, 
                                        Matrix(double) *sp_tr, 
                                        Matrix(double) *n,
                                        double ep_n)
{
  double ksi_nrm = compute_normal(J2P,mat_e,s_tr,sp_tr,n);
  return ksi_nrm - sqrt(2.0/3.0)*(J2P->k0 + (J2P->beta)*(J2P->hp)*ep_n);  
}

int J2_plasticity_integration_alg(double *sp_out,
                                  double *ep_out,
                                  double *gamma_out,
                                  double *F_in, 
                                  double *Fn_in,
                                  double *sp_n_in,
                                  double ep_n,
                                  MATERIAL_J2_PLASTICITY *J2P,
                                  MATERIAL_ELASTICITY *mat_e)
{
  int err = 0;
  
  double G = mat_e->G;

  Matrix(double) F,Fn,sp,sp_n;
     F.m_row =    F.m_col = DIM_3;    F.m_pdata = F_in;  
    Fn.m_row =   Fn.m_col = DIM_3;   Fn.m_pdata = Fn_in;
    sp.m_row =   sp.m_col = DIM_3;   sp.m_pdata = sp_out;  
  sp_n.m_row = sp_n.m_col = DIM_3; sp_n.m_pdata = sp_n_in;  

   
    
  Matrix(double) bbar, s_tr, sp_tr, n, s0;
  Matrix_construct_init(double, bbar,DIM_3,DIM_3,0.0);  
  Matrix_construct_init(double, s_tr,DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,sp_tr,DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,    n,DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,   s0,DIM_3,DIM_3,0.0);  
  
  double J = 0.0;
  Matrix_det(F, J);
  double J23 = pow(J,-2.0/3.0);
    
  // compute bbar at n + 1
  err += compute_bbar_of_J(&bbar, &F,J);
  
  // compute mu_bar
  double mu_bar = G*J23*(Mat_v(bbar,1,1) + Mat_v(bbar,2,2) + Mat_v(bbar,3,3))/3.0;  
      
  // compute sp_tr
  err += compute_sp_tr(&F,&Fn,&sp_n,&sp_tr);
                        
  // compute s_tr
  err += compute_s0(G, &bbar, &s0);
  Matrix_AplusB(s_tr,1.0,s0,-1.0,sp_tr);

  // compute ksi_tr and the normal of plastic loading
      
  double phi = phi_yield_functioncompute_normal(J2P,mat_e,&s_tr,&sp_tr,&n, ep_n);

  if(phi <= J2P_INT_ALG_TOL)
  {
    *gamma_out = 0.0;
    *ep_out = ep_n;
    memcpy(sp_out, sp_tr.m_pdata, DIM_3x3*sizeof(double));
  } 
  else 
  {
    *gamma_out = phi/(2.0* mu_bar*(1.0 + J2P->hp/(3.0*G)*(1.0 - J2P->beta)
                            + J2P->beta*J2P->hp/(3.0*mu_bar)));
    double tmp = 2.0*mu_bar*(*gamma_out);
    *ep_out = ep_n + sqrt(2.0/3.0) *(*gamma_out);
    
    for(int i = 0; i < DIM_3x3; i++) 
      sp_out[i] = sp_tr.m_pdata[i] + tmp*n.m_pdata[i];
  }  
  
  Matrix_cleanup(bbar);
  Matrix_cleanup(s_tr);
  Matrix_cleanup(sp_tr);
  Matrix_cleanup(n);
  Matrix_cleanup(s0);
  return err;	
}

int compute_S0_Sbar(Matrix(double) *S0,
                    Matrix(double) *Sbar,
                    Matrix(double) *F,
                    Matrix(double) *sp,
                    ELASTICITY *elast)
{
  int err = 0;
  // compute the current configuration deviatoric stresses
  
  double G = (elast->mat)->G;
  double kappa = (elast->mat)->kappa;

  Matrix(double) s0,sbar,bbar, C, CI;
  Matrix_construct_init(double,s0,  DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,sbar,DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,bbar,DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,C,   DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,CI,  DIM_3,DIM_3,0.0);    
  
  err += compute_bbar(&bbar, F);
  err += compute_s0(G, &bbar,&s0);  
  Matrix_AplusB(sbar,1.0,s0,-1.0, *sp);

  //perform pull-back */
  err += compute_pull_back(F,&s0,S0);
  err += compute_pull_back(F,&sbar,Sbar);                        

  //compute volumetric stress
  Matrix_AxB(C,1.0,0.0,*F,1,*F,0);
  Matrix_inv(C,CI);
  double J = 0.0;
  Matrix_det(*F,J);
  double dudj = 0.0;
  elast->compute_dudj(&dudj,J);
  for (int i = 0; i < DIM_3x3; i++) 
  {
    double tmp = kappa*J*dudj*CI.m_pdata[i];
    S0->m_pdata[i] += tmp;
    Sbar->m_pdata[i] += tmp;
  }

  Matrix_cleanup(s0);   
  Matrix_cleanup(sbar);
  Matrix_cleanup(bbar);
  Matrix_cleanup(C);
  Matrix_cleanup(CI);  
  return err;
}

// compute the deviatoric initial/unloading tangent in the reference
// configuration
int compute_unloading_Aep_dev(Matrix(double) *Aep_dev,
                              Matrix(double) *F,
                              Matrix(double) *Fn,
                              Matrix(double) *sp_n,
                              double G)
{
  int err = 0;
  Matrix_init(*Aep_dev, 0.0);

  Matrix(double) C,CI,Spn,eye;
  Matrix_construct_init(double,C,   DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,CI,  DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,Spn, DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,eye, DIM_3,DIM_3,0.0);
  
  Matrix_eye(eye,3);
  
  // compute C and related terms
  //compute volumetric stress
  Matrix_AxB(C,1.0,0.0,*F,1,*F,0);
  Matrix_inv(C,CI);
  double J23 = 0.0;
  Matrix_det(C,J23);
  J23 = pow(J23, -1.0/3.0);
  
  double Cpp = Mat_v(C,1,1) + Mat_v(C,2,2) + Mat_v(C,3,3);

  // compute pull-back of spn
  err += compute_pull_back(F,sp_n,&Spn);

  double CSp = 0.0;
  for(int i = 0; i < DIM_3x3; i++) 
    CSp += C.m_pdata[i]*Spn.m_pdata[i];

  for (int i = 1; i <= DIM_3; i++)
  {
    for (int j = 1; j <= DIM_3; j++)
    {
      for (int k = 1; k <= DIM_3; k++) 
      {
        for (int l = 1; l <= DIM_3; l++)
        {
          Tns4_v(*Aep_dev,i,j,k,l) = ((2.0/3.0*J23*(G*Cpp-CSp)*(Mat_v(CI,i,k)*Mat_v(CI,j,l) 
                                     + Mat_v(CI,i,j)*Mat_v(CI,k,l)/3.0))
                                     - 2.0/3.0*J23*(Mat_v(CI,i,j)*(G*Mat_v(eye,k,l) - Mat_v(Spn,k,l))
                                           + (G * Mat_v(eye,i,j) - Mat_v(Spn,i,j))*Mat_v(CI,k,l)));
        }
      }
    }
  }

  Matrix_cleanup(C);
  Matrix_cleanup(CI);
  Matrix_cleanup(Spn);  
  Matrix_cleanup(eye);  
  return err;
}


// compute the deviatoric plastic loading tangent in the reference
// configuration
int compute_loading_Aep_dev(Matrix(double) *Aep_dev,
                            Matrix(double) *F,
                            Matrix(double) *Fn,
                            Matrix(double) *sp_n,
                            double gamma,
                            MATERIAL_J2_PLASTICITY *J2P,
                            MATERIAL_ELASTICITY *mat_e)
{
  int err = 0;
  double G = mat_e->G;

  Matrix_init(*Aep_dev, 0.0);

  Matrix(double) C,CI,Spn,eye,bbar,devbbar,FI,sp_tr,Fubar,s_tr,normal,normal2,zeros,aep;
  Matrix_construct_init(double,C,       DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,CI,      DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,Spn,     DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,eye,     DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,bbar,    DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,devbbar, DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,FI,      DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,sp_tr,   DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,Fubar,   DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,s_tr,    DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,normal,  DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,normal2, DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,zeros,   DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,aep,     DIM_3x3x3x3,1,0.0);
  
  Matrix_eye(eye,3);
  
  Matrix_AxB(C,1.0,0.0,*F,1,*F,0);
  Matrix_inv(C,CI);
  double J23 = 0.0;
  Matrix_det(*F,J23);
  J23 = pow(J23, -2.0/3.0);
  
  Matrix_inv(*F,FI);
  err += compute_bbar(&bbar, F);
  err += compute_dev(&bbar, &devbbar);  
  Matrix_inv(*F,FI);

  // compute Itr = G tr(bbar) - tr(Fubar spn Fubar')
  double Itr = 0;
  {
    err += compute_Fubar(F,Fn,&Fubar);
    err += compute_push_forward(&Fubar,sp_n,&sp_tr);    
    double trace_sp = sp_tr.m_pdata[0] + sp_tr.m_pdata[4] + sp_tr.m_pdata[8];
    sp_tr.m_pdata[0] -= trace_sp/3.0;
    sp_tr.m_pdata[4] -= trace_sp/3.0;
    sp_tr.m_pdata[8] -= trace_sp/3.0;
    Itr = G*(bbar.m_pdata[0] + bbar.m_pdata[4] + bbar.m_pdata[8]) - trace_sp;
  }

  // compute s_tr, normal and ||s_tr||
  err += compute_s0(G, &bbar,&s_tr);
  double norm_s_tr = 0.0;
  for (int i = 0; i < DIM_3x3; i++) {
    s_tr.m_pdata[i] -= sp_tr.m_pdata[i];
    norm_s_tr += s_tr.m_pdata[i]*s_tr.m_pdata[i];
  }
  norm_s_tr = sqrt(norm_s_tr);

  double ksi_nrm = compute_normal(J2P,mat_e,&s_tr,&sp_tr,&normal);
  Matrix_AxB(normal2,1.0,0.0,normal,0,normal,0);
                      
  // compute factors
  double mu_bar = G*J23*(bbar.m_pdata[0] + bbar.m_pdata[4] + bbar.m_pdata[8])/3.0;
  double f0 = 1.0 - 2.0*mu_bar*gamma/norm_s_tr;
  double del0 = 1.0 + J2P->hp/(3.0*mu_bar);
  double f1 = (1.0/del0 - 1.0 + f0);
  double del1 = f1 * 2.0/3.0*Itr;
  double del2 = (1.0/del0 - 1.0)*4.0/3.0*mu_bar*gamma - 2.0/3.0*norm_s_tr*f1;
  double del3 = 2.0*norm_s_tr*f1;
  double del4 = (1.0/del0 - 1.0)*4.0/3.0*G*gamma*J23;

  // compute aep according to box 5 in Simo and Ju 1989
  for (int i = 1; i <= DIM_3; i++)
  {
    for (int j = 1; j <= DIM_3; j++)
    {
      for (int k = 1; k <= DIM_3; k++)
      {
        for (int l = 1; l <= DIM_3; l++)
        {
          Tns4_v(aep,i,j,k,l) = (f0*(2.0/3.0*Itr*(Mat_v(eye,i,k)*Mat_v(eye,l,j) - Mat_v(eye,i,j)*Mat_v(eye,k,l)/3.0)
                                - 2.0/3.0*(Mat_v(s_tr,i,j)*Mat_v(eye,k,l) + Mat_v(eye,i,j)*Mat_v(s_tr,k,l)))
                       - del1*Mat_v(normal,i,j)*Mat_v(normal,k,l)
                       - del2*(Mat_v(normal,i,j)*Mat_v(eye,k,l) + Mat_v(normal,k,l)*Mat_v(eye,i,j))*0.5
                       - del3*(Mat_v(normal,i,j)*Mat_v(normal2,k,l) + Mat_v(normal,k,l)*Mat_v(normal2,i,j))*0.5
                       - del4*(Mat_v(normal,i,j)*Mat_v(devbbar,k,l) + Mat_v(normal,k,l)*Mat_v(devbbar,i,j))*0.5);
        }
      }
    }
  }
  
  err += compute_pull_Tensor4(&FI,&aep,Aep_dev);

  Matrix_cleanup(C);
  Matrix_cleanup(CI);
  Matrix_cleanup(Spn);
  Matrix_cleanup(eye);
  Matrix_cleanup(bbar);
  Matrix_cleanup(devbbar);
  Matrix_cleanup(FI);
  Matrix_cleanup(sp_tr);
  Matrix_cleanup(Fubar);
  Matrix_cleanup(s_tr);
  Matrix_cleanup(normal);
  Matrix_cleanup(normal2);
  Matrix_cleanup(zeros);
  Matrix_cleanup(aep);
    
  return err;
}

int compute_Lbar(Matrix(double) *Lbar,
                 Matrix(double) *F,
                 Matrix(double) *Fn,
                 Matrix(double) *sp_n,
                 double gamma,
                 MATERIAL_J2_PLASTICITY *J2P,             
                 ELASTICITY *elast)
{
  int err = 0;
  double kappa = (elast->mat)->kappa;
  double G = (elast->mat)->G;
  double J = 0.0;
  Matrix_det(*F,J);
  
  Matrix(double) C, CI;
  Matrix_construct_init(double,C,  DIM_3,DIM_3,0.0);
  Matrix_construct_init(double,CI, DIM_3,DIM_3,0.0); 
  
  Matrix_AxB(C,1.0,0.0,*F,1,*F,0);
  Matrix_inv(C,CI);
    
  double dudj = 0.0;
  double d2udj2 = 0.0;

  /* compute deviatoric tangent */
  if (gamma > 0)
    err += compute_loading_Aep_dev(Lbar,F,Fn,sp_n,gamma,J2P,elast->mat);
  else 
    err += compute_unloading_Aep_dev(Lbar,F,Fn,sp_n,G);    

  elast->compute_dudj(&dudj,J);
  elast->compute_d2udj2(&d2udj2,J);

  double coeff_1 = kappa*J*(dudj+J*d2udj2);
  double coeff_2 = 2.0*kappa*J*dudj;
    
  for (int i=1; i <= DIM_3; i++) 
  {
    for (int j=1; j <= DIM_3; j++) 
    {
      for (int k=1; k <= DIM_3; k++) 
      {
	      for (int l=1; l <= DIM_3; l++) 
	      {
	        // Deviatoric + Volumetric stiffness
	        Tns4_v(*Lbar,i,j,k,l) += ((coeff_1*Mat_v(CI,i,j)*Mat_v(CI,k,l))
                                -(coeff_2*Mat_v(CI,i,k)*Mat_v(CI,l,j)));
	      }
      }
    }
  }
  
  Matrix_cleanup(C);
  Matrix_cleanup(CI); 
  return err;
}

int J2_plasticity_update_elasticity(MATERIAL_J2_PLASTICITY *J2P,
                                    ELASTICITY *elast,
                                    double *F_in,
                                    double *Fn_in,
                                    double *sp_in,
                                    double *sp_n_in,
                                    double gamma,
                                    const int compute_stiffness)
{
  int err = 0;
  Matrix(double) F,Fn,sp,sp_n,S0,S,L;
     F.m_row =    F.m_col = DIM_3;    F.m_pdata = F_in;
    Fn.m_row =   Fn.m_col = DIM_3;   Fn.m_pdata = Fn_in;
    sp.m_row =   sp.m_col = DIM_3;   sp.m_pdata = sp_in;    
  sp_n.m_row = sp_n.m_col = DIM_3; sp_n.m_pdata = sp_n_in;
     S.m_row =    S.m_col = DIM_3;    S.m_pdata = elast->S;
     L.m_row = DIM_3x3x3x3;
     L.m_col = 1;  
     L.m_pdata = elast->L;
         
  Matrix_construct_init(double,S0,DIM_3,DIM_3,0.0);
  err += compute_S0_Sbar(&S0,&S,&F,&sp,elast);
  
  if(compute_stiffness)
    err += compute_Lbar(&L,&F,&Fn,&sp_n,gamma,J2P,elast); //compute stiffness

  Matrix_cleanup(S0);
  return err;
}