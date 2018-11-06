#include "test_functions.h"

double **product_matrix_matrix(double **matrixA, double **matrixB, int num) {
  int i, j, k;
  double  **new_matrix = allocate_matrix(num, num);
  zero_matrix(new_matrix, num, num);

  for (i = 0; i < num; i++){
    for (j = 0; j < num; j++){
      for (k = 0; k < num; k++){
	new_matrix[i][j] = new_matrix[i][j] + matrixA[i][k] * matrixB[k][j];
      }
    }
  }

  return new_matrix;
}
