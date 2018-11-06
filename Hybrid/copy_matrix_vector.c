#include "test_functions.h"

void copy_matrix_vector(double **matrix, double *vector, int num, int col, int num_thrds) {
  int  j;
  zero_vector(vector, num, num_thrds);
  #pragma omp parallel for private(j) num_threads(num_thrds)
  for (j = 0; j < num; j++){
    vector[j] = matrix[col][j+1];
  }
}
