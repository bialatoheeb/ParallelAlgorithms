#include<stdlib.h>
#include <test_functions.h>

double *product_vector_scalar(double *vector, int num, double mult){
  int i;
  double *new_vector = allocate_vector(num);

  for (i=0; i< num; i++){
    new_vector[i] = mult * vector[i];
  }
  return new_vector;
}

