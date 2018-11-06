#include "test_functions.h"

void print_matrix(double **A, int num1, int num2){

  int i,j;
  for(i =0; i < num1; i++){
    for(j=0; j < num2; j++){
      printf("%.3e\t", A[i][j]);
    }
    printf("\n");
  }
  printf("\n\n");

}
