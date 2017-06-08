#include <iostream>
#include <cmath>
#include "data_structure.h"

using namespace std;

#define TOL 1.0e-15

/// print test statement
///
/// \param[in] test_no test number
/// \param[in] str printing statement
/// \return non-zero on internal error
int print_test_header(const int test_no, const char *str)
{
  int err = 0;
  cout << "==================================================\n";
  cout << "Test" <<  test_no << ": " << str;
  cout << "--------------------------------------------------\n";
  return err;
}

/// print data array
///
/// \param[in] A Matrix<double> object
/// \return non-zero on internal error
int print_data(Matrix<double> &A, const char *name)
{
  int err = 0;
  cout << name << " = [\n";
  for(int i=1; i<=A.m_row; i++)
	{
		for(int j=1; j<=A.m_col; j++)
		  cout << A.get_a_value(i, j) << " ";

		if(i==A.m_row)
		  cout << "];\n";
		else
		  cout << "\n";
	}
	return err;
}

/// determine whether test is good or failed
///
/// \param[in] A Matrix<double> object
/// \param[in] data double array
/// \return non-zero on internal error
int judge_test(Matrix<double> &A, const double *data)
{
  int err = 0;
  double diff = 0.0;
  if((A.m_row*A.m_col == 0) || (A.m_pdata == NULL))
  {  
    cout << "WARNING!, Test failed\n\n";
    err++;    
    return err;
  }
    
  for(int ia = 0; ia<A.m_row*A.m_col; ia++)
    diff = (A.m_pdata[ia]-data[ia])*(A.m_pdata[ia]-data[ia]);
  
  double norm = sqrt(diff);
  if(norm<TOL)
    cout << ">> Test passed\n\n";
  else
    cout << ">> WARNING!, Test faild\n\n";

  return err;  
}

int basic_data_management(int test_no)
{
  int err = 0;  
  double data[9] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};

  // test : Construct data cell with a data pointer
  err += print_test_header(test_no++, "Construct a matrix with a data array\n");
  Matrix<double> A(3,3,data);
  err += print_data(A, "A");
  err += judge_test(A,data);
  
  // test : B = A
  err += print_test_header(test_no++, "Perform A = B\n");
  Matrix<double> B;
  B = A;
  err += print_data(B, "B");
  err += judge_test(B,data);
  
  // test : C[4 by 4] = A[3 by 3]
  err += print_test_header(test_no++, "Perform C[4 by 4] = A[3 by 3], re-size matrix\n");
  Matrix<double> C(4,4,1.0);
  
  err += print_data(C, "C");
  err += print_data(A, "A");
  C = A;
  err += print_data(C, "C");
  err += judge_test(C,data);
    
  return test_no;
} 

int matrix_operation(int test_no)
{
  int err = 0;
  double data[9] = {0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3};
  
  // test : test for A = B + C
  err += print_test_header(test_no++, "A = B + C\n");
  Matrix<double> A, B, C;
  
  B.initialization(3,3, 0.1);
  C.initialization(3,3, 0.2);
  
  A = B + C;
  A.print("A");
  B.print("B");
  C.print("C");
  err += judge_test(A,data);
  
  // test : test for A = B - C
  B.initialization(3,3, 0.5);
  err += print_test_header(test_no++, "A = B - C\n");
    
  A = B - C;
  A.print("A");
  B.print("B");
  C.print("C");
  err += judge_test(A,data);
  
  // test : test for A = B*C
  B.initialization(3,2, 0.5);
  C.initialization(2,3, 0.3);  
  err += print_test_header(test_no++, "A = B * C\n");
    
  A = B*C;
  A.print("A");
  B.print("B");
  C.print("C");
  err += judge_test(A,data);
  
  // test : test for B *= C
  err += print_test_header(test_no++, "B *= C\n");

  B.print("B");
  C.print("C");  
  B.prod(C);
  B.print("B");

  err += judge_test(B,data);
  
  // test : test for A = a*B + C;
  err += print_test_header(test_no++, "A = a*B + C\n");
  
  double a = 0.2;
  B.initialization(3,3,1.0);
  C.initialization(3,3,0.1);
 
  A = a*B + C;

  printf("a = %f\n", a);
  B.print("B");
  C.print("C");

  A.print("A");

  err += judge_test(A,data);   
  
  return test_no;
} 

int main(int argc,char *argv[])
{
  int err = 0;  
  int test_no = 0;

  test_no += basic_data_management(test_no);
  test_no += matrix_operation(test_no);
  return err;
}
