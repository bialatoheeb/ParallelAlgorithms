#include "test_functions.h"

void print_vector(double *A, int num1){

  int i;
  for(i =0; i < num1; i++){
    printf("%.3e\t", A[i]);
  }
  printf("\n\n");
  
}
