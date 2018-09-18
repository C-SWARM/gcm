#ifndef __H__GCM_DATA_STRUCTURE_H__
#define __H__GCM_DATA_STRUCTURE_H__

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "mkl_cblas.h"
#include "math_help.h"
#include <math.h>

using namespace std;

namespace gcm {
#define DEBUG_PRINT_GCM_DATA_STRUCTURE 0

#define CRM CblasRowMajor
#define CNT CblasNoTrans
#define CYT CblasTrans

/// C = aAxB + bC
/// C[m,n] = a*A[m,k] x B[k,n] + b*C[m,n]
// cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,a,A,k,B,n,b,C,n);
template <class T>
int Matrix_AxB_large(const T &C, 
                     double a, 
                     double b, 
                     const T &A, 
                     const int AT, 
                     const T &B, 
                     const int BT)
{
  int err = 0;               
  if(AT==0 && BT==0)
    cblas_dgemm(CRM,CNT,CNT,A.m_row,B.m_col,A.m_col,a,A.m_pdata,A.m_col,B.m_pdata,B.m_col,b,C.m_pdata,C.m_col);
    
  if(AT==1 && BT==0)
    cblas_dgemm(CRM,CYT,CNT,C.m_row,C.m_col,B.m_row,a,A.m_pdata,A.m_col,B.m_pdata,B.m_col,b,C.m_pdata,C.m_col);

  if(AT==0 && BT==1)
    cblas_dgemm(CRM,CNT,CYT,C.m_row,C.m_col,A.m_col,a,A.m_pdata,A.m_col,B.m_pdata,B.m_col,b,C.m_pdata,C.m_col);

  if(AT==1 && BT==1)
    cblas_dgemm(CRM,CYT,CYT,C.m_row,C.m_col,A.m_row,a,A.m_pdata,A.m_col,B.m_pdata,B.m_col,b,C.m_pdata,C.m_col);

  return err;
}

/// check Matrix size A and B for operatations: multiplication
///
/// if not same, print function name and error massage.
///
/// \param[in] A 1st matrix to be compared
/// \param[in] B 2nd matrix to be compared
/// \param[in] func_name function name calling this function
/// \param[in] lineno line number
/// \return true or false
template <class T>
bool check_matrix_size_AxB(const T &A, const T &B, const char *func_name, const int lineno)
{
  bool is_it_same = true;
  if(A.m_col !=B.m_row)
  {  
    cout << "Matrix<T>" << func_name << "(" << lineno <<"): Error Matrix size does not match. ";
    cout << "[" << A.m_row << "x" << A.m_col <<"] * ";
    cout << "[" << B.m_row << "x" << B.m_col <<"] \n";    
    is_it_same = false;
  }
  return is_it_same;
}

/// check Matrix size A, B and C for operatations: multiplication A*B*C
///
/// if not same, print function name and error massage.
///
/// \param[in] A 1st matrix to be compared
/// \param[in] B 2nd matrix to be compared
/// \param[in] C 3rd matrix to be compared
/// \param[in] func_name function name calling this function
/// \param[in] lineno line number
/// \return true or false
template <class T>
bool check_matrix_size_AxBxC(const T &A, const T &B, const T &C, const char *func_name, const int lineno)
{
  bool is_it_same = true;
  if((A.m_col !=B.m_row) || (B.m_col != C.m_row))
  {  
    cout << "Matrix<T>" << func_name << "(" << lineno <<"): Error Matrix size does not match. ";
    cout << "[" << A.m_row << "x" << A.m_col <<"] * ";
    cout << "[" << B.m_row << "x" << B.m_col <<"] * ";
    cout << "[" << C.m_row << "x" << C.m_col <<"] \n";        
    is_it_same = false;
  }
  return is_it_same;
}

/// check Matrix size A and B for operatations: =, +, -
///
/// if not same, print function name and error massage.
///
/// \param[in] A 1st matrix to be compared
/// \param[in] B 2nd matrix to be compared
/// \param[in] func_name function name calling this function
/// \param[in] lineno line number
/// \return true or false
template <class T>
bool check_matrix_size_A_B(const T &A, const T &B, const char *func_name, const int lineno)
{
  bool is_it_same = true;
  if((A.m_row !=B.m_row)||(A.m_col !=B.m_col))
  {  
    cout << "Matrix<T>" << func_name << "(" << lineno <<"): Error Matrix size does not match!\n";
    is_it_same = false;
  }
  return is_it_same;
}

/// check Matrix size A and B for operatations: =, +, -
///
/// if not same, print function name and error massage.
///
/// \param[in] A 1st matrix to be compared
/// \param[in] B 2nd matrix to be compared
/// \param[in] func_name function name calling this function
/// \param[in] lineno line number
/// \return true or false
template <class T>
bool check_square_matrix(const T &A, const char *func_name, const int lineno)
{
  bool is_it_square = true;
  if((A.m_row == 0) || (A.m_row !=A.m_row))
  {
    cout << "Matrix<T>" << func_name << "(" << lineno <<"): Error Matrix is not square!\n";
    is_it_square = false;
  }
  return is_it_square;
}

/// GCM_DATA_STRUCTURE class 
///
/// Topest class of the GSM datastructure
/// Matrix type of data array for any data types 
template <class T> class Matrix
{
public:
  
  /// number of row
  int m_row;
  /// number of col 
  int m_col;
  /// memory for data array
  T *m_pdata;
  /// check for reference use 
  bool is_matrix_created;

  // constructors
  Matrix();                                
  Matrix(const Matrix<T> &A);
  Matrix(const int m_size, const int n_size);
  Matrix(const int m_size, const int n_size, const T &value);
  Matrix(const int m_size, const int n_size, const T *value);
  void set_nulls(void);

  // destructor
  ~Matrix();
  
  // initialize Matrix
  void initialization(const int m_size, const int n_size, const T &value);       
  void initialization(const int m_size, const int n_size);                       
  void initialization(const Matrix<T> &A);                           
  void initialization(const int m_size, const int n_size, const T *p);
  
  void use_reference(const int m_size, const int n_size, T *p)
  {
    m_row = m_size;
    m_col = n_size;
    m_pdata = p;
  }

  // deallocate memory of member array  
  void cleanup(void);
  
  // check data array size
  int check_size_and_null(const int m, const int n);             

  // set a value 
  void set_a_value(const int m, const int n, const T &value);      
  void set_values(const int m_size, const int n_size, const T *p);  
  void set_values(const T &p);

  // get a value
  T  get_a_value(const int m, const int n);
  T* get_a_value_pointer(const int m, const int n);
  
  /// indexing data through bracket
  ///
  /// e.g. A(a, b) = c;
  ///
  /// \param[in] a row index
  /// \param[in] b column index
  /// \return data component at a, b
  constexpr T& operator () (const int a, const int b)
  {
    return *(m_pdata + a*m_col + b);
  }
  

  /// indexing data through bracket
  ///
  /// e.g. A(a) = c;
  ///
  /// \param[in] a row index and assumed that number of column 1
  /// \return data component at a, b
  constexpr T& operator () (const int a)
  {
    return *(m_pdata+a);
  }
  
  /// operator overriding for Matrix  
  ///
  /// \param[in] obj a matrix to be assigned
  /// \return this matrix
  Matrix<T> & operator = (const Matrix<T> &obj)
  {
    initialization(obj);
    return *this;
  }

  // add  
  void add(const Matrix<T> &B);						          // this += B
  void add(const Matrix<T> &A, const Matrix<T> &B); // this = A + B
  void add(const double b);                         // this += b  
  void add(const Matrix<T> &A, const double b);     // this = A + b
  
  // subtract
  void sub(const Matrix<T> &B);	                    // this -= B
  void sub(const Matrix<T> &A, const Matrix<T> &B);	// this = A - B
  void sub(const double b);	                        // this = b - this 
  void sub(const double a, const Matrix<T> &B);	    // this = a - B

  // product
  void prod(const Matrix<T> &B);                     // this *= B
  void prod(const Matrix<T> &A, const Matrix<T> &B); // this = A*B
  void prod(const Matrix<T> &A, const Matrix<T> &B, const Matrix<T> &C); // this = A*B*C
  void prod(const Matrix<T> &A, const double b);     // this = A*b
  void prod(const double b);                         // this *= b (scalar)
  void prod(const Matrix<T> &A, 
            const int AT, 
            const Matrix<T> &B,
            const int BT);
  
  void trans(void);                      // transpose
  void trans(Matrix<T> &A);
  double norm(void);                     // sqrt(this:this)
  double ddot(Matrix<T> &B);             // this:B
  double ox(Matrix<T> &A, Matrix<T> &B); // this = A /otimes B 
  void input_from_str(char *str);
  void print(void);
  void print(const char *name, const char *format = "%e");
  int inv(void);
  int inv(const Matrix<T> &A);
 };

/// default constructor
template <class T> Matrix<T>:: Matrix()
{
  set_nulls();  
}

/// constructor with a already exist object
///
/// \param[in] A a Matrix
/// \return no return
template <class T> Matrix<T>:: Matrix(const Matrix<T> &A)
{
  set_nulls();
  this->initialization(A);
}

/// constructor with initialization of m by n matrix
///
/// \param[in] m_size row index
/// \param[in] n_size column index
/// \return no return
template <class T> Matrix<T>:: Matrix(const int m_size, 
                                     const int n_size)
{
  set_nulls();
  this->initialization(m_size, n_size);
}

/// constructor with initialization of m by n Matrix and fill data with a (type T) value
///
/// \param[in] m_size row index
/// \param[in] n_size column index
/// \param[in] value a data (type T) value
/// \return no return
template <class T> Matrix<T>:: Matrix(const int m_size, 
                                     const int n_size, 
                                     const T &value)
{
  set_nulls();
  this->initialization(m_size, n_size, value);
}

/// constructor with initialization of m by n Matrix and fill data with a (type T) data array
///
/// \param[in] m_size row index
/// \param[in] n_size column index
/// \param[in] value a data (type T) value
/// \return no return
template <class T> Matrix<T>:: Matrix(const int m_size, 
                                     const int n_size, 
                                     const T *value)
{
  set_nulls();
  this->initialization(m_size, n_size, value);
}

/// default destructor
template <class T> Matrix<T>:: ~Matrix()
{
  cleanup();
}

/// set all members to NULL
///
/// \return void
template <class T> void Matrix<T>::set_nulls(void)
{
  m_row=0;
  m_col=0;
  m_pdata=NULL;
  is_matrix_created = false;
  if(DEBUG_PRINT_GCM_DATA_STRUCTURE)
    cout << "Debug mode >> set NULL ...\n";
}

/// initialize Matrix to have m by n data array and fill the data array with (type T) value
///
/// \param[in] m_size row index
/// \param[in] n_size column index
/// \param[in] value a data (type T) value
/// \return void
template <class T> void Matrix<T>::initialization(const int m_size, 
                                                  const int n_size, 
                                                  const T &value)
{
  this->initialization(m_size, n_size);
  
  for(int i=0; i<m_row; i++)
  {
    for(int j=0; j<m_col; j++)
    {
      set_a_value(i, j, value);
    }
  }
}

/// initialize Matrix to have m by n data array
///
/// \param[in] m_size row index
/// \param[in] n_size column index
/// \return void
template <class T> void Matrix<T>::initialization(const int m_size, 
                                                  const int n_size)
{
  // if size of data array is same as inputs
  // no need to do memory allocation.

  if(DEBUG_PRINT_GCM_DATA_STRUCTURE)
    cout << "Debug mode >> construct [" << m_size << "x" << n_size << "] data cell...\n";

  if((m_row*m_col == m_size*n_size) && (m_pdata != NULL)) 
  {
    if(DEBUG_PRINT_GCM_DATA_STRUCTURE)
      cout << "Debug mode >> already constructed.done.\n";

    m_row = m_size; // just set the size     
    m_col = n_size;    
    return;
  }

  if(m_pdata != NULL)
  {
    if(DEBUG_PRINT_GCM_DATA_STRUCTURE)
    {  
      cout << "Debug mode >> already constructed but size is different: ";
      cout << "Debug mode >> this  = [" << m_row  << "x" << m_col  << "] vs "; 
      cout << "Debug mode >> input = [" << m_size << "x" << n_size << "], "; 
      cout << "cleanup...\n";
    }
    cleanup();
  }
    
  m_row = m_size;
  m_col = n_size;
  m_pdata = new T [m_size*n_size] ();
  is_matrix_created = true;
    
  if(m_pdata == NULL)
    cout << "Memory is not allocated. " << __func__ << ":" << __FILE__ << ":" << __LINE__ << "\n";
  else if(DEBUG_PRINT_GCM_DATA_STRUCTURE)
    cout << "Debug mode >> [" << m_row << "x" << m_col << "] data cell is constructed.done.\n";    
}

/// initialize Matrix with a already exist object
///
/// If the size of this Matrix is different from the input object,
/// this will be re-sized.
///
/// \param[in] A a Matrix object
/// \return void
template <class T> void Matrix<T>::initialization(const Matrix<T> &A)
{
  this->initialization(A.m_row, A.m_col);
  for(int ia=0; ia<A.m_row*A.m_col; ia++)
    m_pdata[ia] = A.m_pdata[ia];     
}

/// initialize Matrix with a data array
///
/// If the size of this Matrix is different than inputs,
/// this will be re-sized.
///
/// \param[in] m_size number of row of input array
/// \param[in] n_size number of columns of input array
/// \param[in] p a data (type T) array to be set
/// \return void
template <class T> void Matrix<T>::initialization(const int m_size, 
                                                  const int n_size, 
                                                  const T *p)
{ 
  this->initialization(m_size, n_size);

  for(int i=0; i<m_row; i++)
  {
    for(int j=0; j<m_col; j++)
      set_a_value(i, j, p[i*m_size + j]);
  }
}

/// deallocate memory and set size to zeros 
///
/// \return void
template <class T> void Matrix<T>::cleanup(void)
{
  if(!is_matrix_created)
    m_pdata = NULL; // release pointers if reference is used

  if(m_pdata!=NULL)
  {
    if(DEBUG_PRINT_GCM_DATA_STRUCTURE)
      cout << "Debug mode >> delete data and set NULL\n";      
    delete[] m_pdata;
    m_pdata = NULL;
    is_matrix_created = false;
  }
  
  m_row = 0;
  m_col = 0;
}

/// check input indices are bounded by the number of row and column
///
/// If indices are outbounded, return 0 otherwise return 1.
///
/// \param[in] m row index
/// \param[in] n column index
/// \return 0: m or n is out of range of the data array
///         1: within the range of the data array
template <class T> int Matrix<T>::check_size_and_null(const int m, 
                                                      const int n)
{
  if((0<m && m<m_row) && (0<n && n<m_col) && (m_pdata != NULL)) 
    return 1;
  else
    return 0;
}

/// set a component(m, n) with a data (type T) value 
///
/// If the size of this Matrix is different from inputs,
/// this will be re-sized.
///
/// \param[in] m row index
/// \param[in] n column index
/// \param[in] value a data (type T) value
/// \return void
template <class T> void Matrix<T>::set_a_value(const int m, 
                                               const int n, 
                                               const T &value)
{
  *(m_pdata + m*m_col + n) = value;
}

/// set all components by a data array
///
/// If the size of this Matrix is different than inputs,
/// this will be re-sized.
///
/// \param[in] m_size number of row of input array
/// \param[in] n_size number of columns of input array
/// \param[in] p a data (type T) array to be set
/// \return void
template <class T> void Matrix<T>::set_values(const int m_size, 
                                              const int n_size, 
                                              const T *p)
{
  initialization(m_size, n_size, p);
}

/// set all components by a data
///
/// If the size of this Matrix is different than inputs,
/// this will be re-sized.
///
/// \param[in] d a data (type T)
/// \return void
template <class T> void Matrix<T>::set_values(const T &d)
{
  this->initialization(m_row, m_col, d);
}

/// get a component(m, n) value
///
/// \param[in] m row index
/// \param[in] n column index
/// \return a data (type T) value
template <class T> T Matrix<T>::get_a_value(const int m, 
                                            const int n)
{
  return *(m_pdata + (m*m_col + n));
}

/// get pointer of a component(m, n)
///
/// If the pointer is used after, this object is destructed,
/// memory assess error could occor.
///
/// \param[in] m row index
/// \param[in] n column index
/// \return a pointer of a data (type T) at m, n
template <class T> T* Matrix<T>::get_a_value_pointer(const int m, 
                                                     const int n)
{
  return m_pdata + (m*m_col + n);
}

/// matrix addition to self
///
/// data calling this function will be updated as a result of
/// operation : this += B
///
/// \param[in] B input matrix to be added
template <class T> void Matrix<T>::add(const Matrix<T> &B)
{
  if((this->m_row==0)||(this->m_col==0))
  {
    this->set_values(B.m_row, B.m_col, B.m_pdata);
    return;
  }
  if((B.m_row==0)||(B.m_col==0))
    return;
    
  if(check_matrix_size_A_B(*this, B, __func__, __LINE__))
  {
    for(int i=0; i<(this->m_row)*(this->m_col); i++)
      this->m_pdata[i] += B.m_pdata[i];

  }
}

/// matrix addition
///
/// data calling this function will be updated as a result of
/// operation : this = A + B
///
/// \param[in] A 1st input matrix to be added
/// \param[in] B 2nd input matrix to be added
template <class T> void Matrix<T>::add(const Matrix<T> &A, const Matrix<T> &B)
{
  if((A.m_row==0)||(A.m_col==0))
    return;
    
  if(check_matrix_size_A_B(A, B, __func__, __LINE__))
  {
    this->initialization(B.m_row, B.m_col);
    for(int i=0; i<(this->m_row)*(this->m_col); i++)
      this->m_pdata[i] = A.m_pdata[i] + B.m_pdata[i];
  }
}

/// matrix addition
///
/// data calling this function will be updated as a result of
/// operation : this += b, b = scalar
///
/// \param[in] b a scalar to be added
template <class T> void Matrix<T>::add(const double b)
{
  for(int i=0; i<(this->m_row)*(this->m_col); i++)
      this->m_pdata[i] += b;
}

/// matrix addition
///
/// data calling this function will be updated as a result of
/// operation : this = A + b, b = scalar
///
/// \param[in] A input matrix to be added
/// \param[in] b a scalar to be added
template <class T> void Matrix<T>::add(const Matrix<T> &A, const double b)
{
  if((A.m_row==0)||(A.m_col==0))
    return;
  
  this->initialization(A.m_row, A.m_col);  
  for(int i=0; i<(this->m_row)*(this->m_col); i++)
      this->m_pdata[i] = A.m_pdata[i] + b;
}

/// matrix subtract to self
///
/// data calling this function will be updated as a result of
/// operation : this -= B
///
/// \param[in] B input matrix to be subtracted
template <class T> void Matrix<T>::sub(const Matrix<T> &B)
{
  if((this->m_row==0)||(this->m_col==0))
  {
    this->set_values(B.m_row, B.m_col, B.m_pdata);
    return;
  }
  if((B.m_row==0)||(B.m_col==0))
    return;
  if(check_matrix_size_A_B(*this, B, __func__, __LINE__))
  {
    for(int i=0; i<(this->m_row)*(this->m_col); i++)
      this->m_pdata[i] -= B.m_pdata[i];

  }
}

/// matrix subtraction
///
/// data calling this function will be updated as a result of
/// operation : this = A - B
///
/// \param[in] A 1st input matrix
/// \param[in] B 2nd input matrix to be subtracted
template <class T> void Matrix<T>::sub(const Matrix<T> &A, const Matrix<T> &B)
{
  if((A.m_row==0)||(A.m_col==0))
    return;
    
  if(check_matrix_size_A_B(A, B, __func__, __LINE__))
  {
    this->initialization(B.m_row, B.m_col);
    for(int i=0; i<(this->m_row)*(this->m_col); i++)
      this->m_pdata[i] = A.m_pdata[i] - B.m_pdata[i];
  }
}

/// matrix subtraction
///
/// data calling this function will be updated as a result of
/// operation : this = b - this
///
/// \param[in] b a scalar
template <class T> void Matrix<T>::sub(const double b)
{
  for(int i=0; i<(this->m_row)*(this->m_col); i++)
    this->m_pdata[i] = b - this->m_pdata[i];
}

/// matrix subtraction
///
/// data calling this function will be updated as a result of
/// operation : this = a - B
///
/// \param[in] a a scalar
/// \param[in] B input matrix to be subtracted
template <class T> void Matrix<T>::sub(const double a, const Matrix<T> &B)
{
  if((B.m_row==0)||(B.m_col==0))
    return;

  this->initialization(B.m_row, B.m_col);
  for(int i=0; i<(this->m_row)*(this->m_col); i++)
    this->m_pdata[i] = a - B.m_pdata[i];
}

/// perform matrix product to itself
///
/// data calling this function will be updated as a result of
/// this *= B
///
/// \param[in] B input matrix to be multiplied
template <class T> void Matrix<T>::prod(const Matrix<T> &B)
{
  if(check_matrix_size_AxB(*this, B, __func__, __LINE__))
  {
    Matrix<T> A(*this);
    this->initialization(A.m_row, B.m_col);
    Matrix_AxB_large(*this,1.0,0.0,A,0,B,0);    
  }  
}        

/// perform two matrices product to itself 
///
/// data calling this function will be updated as a result of
/// this = A*B
/// \param[in] A 1st input matrix to be multiplied
/// \param[in] B 2nd input matrix to be multiplied
template <class T> void Matrix<T>::prod(const Matrix<T> &A, const Matrix<T> &B)
{
  if(check_matrix_size_AxB(A, B, __func__, __LINE__))
  {
    this->initialization(A.m_row, B.m_col);
    Matrix_AxB_large(*this,1.0,0.0,A,0,B,0);
  }
  else
    abort(); 
}

/// perform three matrices product to itself 
///
/// data calling this function will be updated as a result of
/// this = A*B*C
/// \param[in] A 1st input matrix to be multiplied
/// \param[in] B 2nd input matrix to be multiplied
/// \param[in] C 3rd input matrix to be multiplied
template <class T> void Matrix<T>::prod(const Matrix<T> &A, const Matrix<T> &B, const Matrix<T> &C)
{
  if(check_matrix_size_AxBxC(A, B, C, __func__, __LINE__))
  {
    this->initialization(A.m_row, C.m_col, (T) 0.0);
    for(int ia=0; ia<A.m_row; ia++)
      for(int id=0; id<C.m_col; id++)
        for(int ib=0; ib<A.m_col; ib++)
          for(int ic=0; ic<B.m_col; ic++)
            this->m_pdata[ia*this->m_col + id] += A.m_pdata[ia*A.m_col + ib]
                                                 *B.m_pdata[ib*B.m_col + ic]
                                                 *C.m_pdata[ic*C.m_col + id];
  }
  else
    abort();
}

/// perform matrices product with a scalar and aply to this 
///
/// data calling this function will be updated as a result of
/// this = A*b
/// \param[in] A input matrix to be multiplied
/// \param[in] b a scalar to be multiplied
template <class T> void Matrix<T>::prod(const Matrix<T> &A, const double b)
{
  this->initialization(A.m_row, A.m_col);
  for(int ia=0; ia<(this->m_row)*(this->m_col); ia++)
    this->m_pdata[ia] = A.m_pdata[ia]*b;
}

/// perform multiply a scalar to itself 
///
/// data calling this function will be updated as a result of
/// this = a*B, a is scalar
/// \param[in] a a scalar value to be multiplied
template <class T> void Matrix<T>::prod(const double value)
{
  for(int ia=0; ia<(this->m_row*this->m_col); ia++)
    this->m_pdata[ia] *= value;    
}

/// perform two matrices product to itself 
///
/// data calling this function will be updated as a result of
/// this = A*B
/// \param[in] A 1st input matrix to be multiplied
/// \param[in] AT if 1, do transpose
/// \param[in] B 2nd input matrix to be multiplied
/// \param[in] BT if 1, do transpose
template <class T> void Matrix<T>::prod(const Matrix<T> &A,
                                        const int AT,
                                        const Matrix<T> &B,
                                        const int BT)
{
  if(check_matrix_size_AxB(A, B, __func__, __LINE__))
  {
    this->initialization(A.m_row, B.m_col);
    Matrix_AxB_large(*this,1.0,0.0,A,AT,B,BT);
  }  
}

/// perform matrix transpose
///
/// data calling this function will be updated as a result of this = A^T
/// data calling this function will be updated as a result of
/// \param[in] A input matrix to be transposed
template <class T> void Matrix<T>::trans(Matrix<T> &A)
{
  if(A.m_row==0 || A.m_col==0)
    return;
    
  if(A.m_row==1 || A.m_col==1)
  {
    this->initialization(A.m_col,A.m_row,A.m_pdata);
    return;
  }
  
  this->initialization(A.m_col,A.m_row);
  for(int ia=0; ia<A.m_row; ia++)
  {
    for(int ib=0; ib<A.m_col; ib++)
      this->set_a_value(ib,ia,A(ia,ib));
  }
}


/// perform matrix transpose
///
/// data calling this function will be updated as a result of this^T
template <class T> void Matrix<T>::trans(void)
{
  if(this->m_row==0 || this->m_col==0)
    return;
        
  if(this->m_row==1 || this->m_col==1)
  {
    int temp = this->m_col;
    this->m_col = this->m_row;
    this->m_row = temp;
    return;
  }
    
  Matrix<T> A(this->m_col,this->m_row);
  for(int ia=0; ia<this->m_row; ia++)
  {
    for(int ib=0; ib<this->m_col; ib++)
      A(ib,ia) = this->get_a_value(ia,ib);
  }
  this->initialization(A);
}

/// perform matrix double contraction
///
/// \param[in] B input matrix to be contracted
/// \return a contracted value of (this:B)
template <class T> double Matrix<T>::ddot(Matrix<T> &B)
{
  if(this->m_row==0 || this->m_col==0)
    return 0.0;
  
  double AB = 0.0;
  if(check_matrix_size_A_B(*this, B, __func__, __LINE__))      
  {
    for(int ia=0; ia<this->m_row*this->m_col; ia++)
      AB += this->m_pdata[ia]*B.m_pdata[ia];
  }
  return AB;
}

/// print Matrix components with input name
///
/// \param[in] name name of matrix
template <class T> void Matrix<T>::print(const char *name, const char *format)
{
  cout << name << " = [\n";

  for(int i=0; i<this->m_row; i++)
  {
    for(int j=0; j<this->m_col; j++)
    {
      printf(format, this->get_a_value(i, j));
      printf(" ");
    }
    cout << "\n";
  }
  cout << "];\n";
}

/// print Matrix components without name input
///
/// default name A will be used
template <class T> void Matrix<T>::print(void)
{
  char A[2] = "A";
  this->print(A);
}

/// compute inverse of Matrix
///
template <class T> int Matrix<T>::inv(void)
{
  if(check_square_matrix(*this, __func__, __LINE__))
  {
    Matrix<T> temp(this->m_row, this->m_col);
    int err = inverse(this->m_pdata, this->m_row, temp.m_pdata);
    memcpy(this->m_pdata, temp.m_pdata, (this->m_row*this->m_col)*sizeof(T));
    return err;

  }
  else
    return 1;
}

/// compute inverse of Matrix
///
template <class T> int Matrix<T>::inv(const Matrix<T> &A)
{
  if(check_square_matrix(A, __func__, __LINE__))
  {
    this->initialization(A.m_row, A.m_col);
    int err = inverse(A.m_pdata, A.m_row, this->m_pdata);
    return err;
  }
  else
    return 1;
}

} // namespace gcm

#endif
