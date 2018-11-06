#include "test_functions.h"

double **add_matrix(double ** matrixA, double ** matrixB, int num1, int num2){
  int i,j;
  double ** new_matrix = allocate_matrix(num1, num2);
  for (i=0; i< num1; i++){
    for (j=0; j < num2; j++){
      new_matrix[i][j] = matrixA[i][j] +  matrixB[i][j];
    }
  }
  return new_matrix;
}

