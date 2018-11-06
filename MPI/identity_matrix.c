#include "test_functions.h"

void identity_matrix(double **matrix, int num) {
  zero_matrix(matrix, num, num);
  int i;
  for (i = 0; i < num; i++)
    matrix[i][i] = 1.0;
}
