#include "test_functions.h"

double *nonlinear_function(double * vector, int num1, int num_thrds){
  double* return_func;
  int i;
  return_func = allocate_vector(num1);
  zero_vector(return_func, num1, num_thrds);
#pragma omp parallel for private(i) num_threads(num_thrds)
  for (i=0; i <num1; i++){
    return_func[i] = pow(vector[i],2);
  }
  return return_func;
}
