#include "test_functions.h"

void permutation_matrix(double **matrix, double ** pivot, int num){
  int i,j,k;
  identity_matrix(pivot, num);
  for (i=0; i < num; i++){
    int max = i;
    for (j = i; j < num; j++){
      if (fabs(matrix[j][i]) > fabs(matrix[max][i]))
	max = j;
    }
    if (max != i){
      for (k=0; k < num; k++){
	int temp = pivot[i][k];
	pivot[i][k] = pivot[max][k];
	pivot[max][k] = temp;

      }
    }

  }
}
