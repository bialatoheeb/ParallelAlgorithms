#include "test_functions.h"

double *add_vector(double *vectorA, double *vectorB, int num, int num_thrds){
  int i;
  double * new_vector = allocate_vector(num);
  
  #pragma omp parallel for private(i) num_threads(num_thrds)
  for (i=0; i< num; i++){
    new_vector[i] = vectorA[i] +  vectorB[i];
  }
  return new_vector;
}
