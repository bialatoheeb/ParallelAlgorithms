#include "test_functions.h"

double **add_matrix(double ** matrixA, double ** matrixB, int num1, int num2, int num_thrds){
  int i,j;
  double ** new_matrix = allocate_matrix(num1, num2);
  
#pragma omp parallel for private(i,j) num_threads(num_thrds)
  for (i=0; i< num1; i++){
    for (j=0; j < num2; j++){
      new_matrix[i][j] = matrixA[i][j] +  matrixB[i][j];
    }
  }
  return new_matrix;
}

