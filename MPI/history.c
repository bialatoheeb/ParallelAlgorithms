#include "test_functions.h"

double*  history(double **Uvec, double **matrixA, int num, double alpha, int start, int end, int col){
  double *history_term;
  double coeff;
  double *temp_vector;
  double *nonlinear_func;
  double * temp1;
  double *temp2;
 
  int i, j;

  history_term = allocate_vector(num);
  zero_vector(history_term, num);
  if (end == 0)
    return history_term;

  temp_vector = allocate_vector(num);
  nonlinear_func = allocate_vector(num);
  temp1 = allocate_vector(num);
  temp2 = allocate_vector(num);
  zero_vector(temp_vector, num);
  zero_vector(nonlinear_func, num);
  zero_vector(temp1, num);
  zero_vector(temp2, num);




  for (i=start; i <= end; i++){
    
    copy_matrix_vector(Uvec, temp_vector, num, i);
    temp1 = product_matrix_vector(matrixA, temp_vector,num);
    nonlinear_func = nonlinear_function(temp_vector, num);
    temp2 = add_vector(temp1, nonlinear_func, num);
    
    if(i == 0)
      coeff = - (col - alpha)*pow(col+1, alpha) + pow(col, alpha)*(2*col - alpha - 1) - pow(col-1, alpha+1);
    else if (i == end && end == col)
      coeff = pow(2, alpha+1) - alpha - 3;
    else
      coeff = pow(col-i+2, alpha+1) - 3*pow(col-i+1, alpha+1) + 3*pow(col-i, alpha+1) - pow(col-i-1, alpha + 1);
    
    temp2 = product_vector_scalar(temp2, num, coeff);
    history_term = add_vector(history_term, temp2, num);
  }
  

  free(temp1);
  free(temp2);
  free(temp_vector);
  free(nonlinear_func);
  return history_term;
}
