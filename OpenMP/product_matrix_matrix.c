#include "test_functions.h"

double **product_matrix_matrix(double **matrixA, double **matrixB, int num, int num_thrds) {
  int i, j, k;
  double  **new_matrix = allocate_matrix(num, num);
  double sum;


#pragma omp parallel for private(i,j,k, sum) num_threads(num_thrds)
  for (i = 0; i < num; i++){
    for (j = 0; j < num; j++){
      sum = 0.0;
      for (k = 0; k < num; k++){
	sum += (matrixA[i][k] * matrixB[k][j]);
      }
      new_matrix[i][j] = sum;
    }
  }

  return new_matrix;
}
