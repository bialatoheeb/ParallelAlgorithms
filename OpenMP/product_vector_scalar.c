#include "test_functions.h"

double *product_vector_scalar(double *vector, int num, double mult, int num_thrds){
  int i;
  double *new_vector = allocate_vector(num);
  #pragma omp parallel for private(i) num_threads(num_thrds)
  for (i=0; i< num; i++){
    new_vector[i] = mult * vector[i];
  }
  return new_vector;
}

