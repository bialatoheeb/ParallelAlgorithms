#include <test_functions.h>


void identity_matrix(double **matrix, int num, int num_thrds) {
  zero_matrix(matrix, num, num, num_thrds);
  int i;
  #pragma omp parallel for private(i) num_threads(num_thrds)
  for (i = 0; i < num; i++)
    matrix[i][i] = 1.0;
}
