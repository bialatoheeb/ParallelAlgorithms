#include "test_functions.h"

double *add_vector(double *vectorA, double *vectorB, int num){
  int i;
  double * new_vector = allocate_vector(num);
  for (i=0; i< num; i++){
    new_vector[i] = vectorA[i] +  vectorB[i];
  }
  return new_vector;
}
