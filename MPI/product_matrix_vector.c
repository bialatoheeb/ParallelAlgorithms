#include "test_functions.h"

double *product_matrix_vector(double **matrixA, double *vector, int num) {
  int i, j, k;
  double  *new_vector = allocate_vector(num);
  double sum = 0;
  zero_vector(new_vector, num);

  for (i = 0; i < num; i++){
    for (j = 0; j < num; j++){
      sum = sum + matrixA[i][j] * vector[j];
    }
    new_vector[i] = sum;
    sum = 0;
  }
  return new_vector;
}
