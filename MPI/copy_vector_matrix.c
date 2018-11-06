#include "test_functions.h"

void copy_vector_matrix(double **matrix, double *vector, int num, int col) {
  int  j;
  for (j = 0; j < num; j++){
    matrix[col][j+1] = vector[j];
  }
}
