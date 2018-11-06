#include "test_functions.h"

void zero_vector(double *vector, int num, int num_threads){
  int j;
  #pragma omp parallele for private(j) num_threads(num_thrds)
  for (j=0; j< num; j++)
    vector[j] = 0.0;
}

