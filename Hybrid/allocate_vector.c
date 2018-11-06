#include "test_functions.h"

double *allocate_vector(int num) {
  double *new_vector = malloc(num *sizeof(double));
  return new_vector;
}
