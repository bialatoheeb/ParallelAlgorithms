#include "test_functions.h"

void copy_matrix_vector(double **matrix, double *vector, int num, int col) {
  int  j;
  zero_vector(vector, num);
  for (j = 0; j < num; j++){
    vector[j] = matrix[col][j+1];
  }
}
