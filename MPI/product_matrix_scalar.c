#include "test_functions.h"

double **product_matrix_scalar(double **matrix, int num1, int num2,  double mult){
  int i,j;
  double **new_matrix = allocate_matrix(num1, num2);

  for (i=0; i< num1; i++){
    for (j=0; j < num2; j++){
      new_matrix[i][j] = mult * matrix[i][j];
    }
  }
  return new_matrix;
}
