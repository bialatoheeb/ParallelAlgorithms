#include "test_functions.h"

double *product_matrix_vector(double **matrixA, double *vector, int num, int num_thrds) {
  int i, j;
  double  *new_vector = allocate_vector(num);
  double sum ;

#pragma omp parallel for private(i, j, sum) num_threads(num_thrds)
  for (i = 0; i < num; i++){
    sum = 0.0;
    for (j = 0; j < num; j++){
      sum += (matrixA[i][j] * vector[j]);
    }
    new_vector[i] = sum;
  }
  
  return new_vector;
}
