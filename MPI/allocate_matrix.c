#include "test_functions.h"

double **allocate_matrix(int num1, int num2) {
  double **new_matrix = (double **) malloc(num1 * sizeof(double*));
  int i;
  for (i = 0; i < num1; i++)
    new_matrix[i] = (double *) malloc(num2 * sizeof(double)) ;
  return new_matrix;
}
