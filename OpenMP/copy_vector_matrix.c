#include "test_functions.h"

void copy_vector_matrix(double **matrix, double *vector, int num, int col, int num_thrds) {
  int  j;
#pragma omp parallel for private(j) num_threads(num_thrds)
  for (j = 0; j < num; j++){
    matrix[col][j+1] = vector[j];
  }
}
