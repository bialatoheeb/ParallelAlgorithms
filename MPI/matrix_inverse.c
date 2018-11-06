#include "test_functions.h"

double ** matrix_inverse(double **matrix, int num){
  double ** L;
  double ** U;
  double ** pivot;
  double ** temp;
  double ** result;
  L = allocate_matrix(num, num);
  U = allocate_matrix(num, num);
  pivot = allocate_matrix(num, num);
  temp = allocate_matrix(num, num);
  result = allocate_matrix(num, num);

  
  LU_decomposition(matrix, L, U, pivot, num);
  temp = forward_substitution(L, pivot, num);
  result = backward_substitution(U, temp, num);


  free_matrix(L, num);
  free_matrix(U, num);
  free_matrix(temp, num);
  free_matrix(pivot, num);
  
  return result;
  
}
