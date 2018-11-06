#include "test_functions.h"

void zero_matrix(double **matrix, int num1, int num2){
  int i, j;
  for (i=0; i< num1; i++){
    for (j=0; j< num2; j++)
      matrix[i][j] = 0.0;
  }
}
