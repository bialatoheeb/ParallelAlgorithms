#include "test_functions.h"

double ** backward_substitution(double ** matrixB, double ** matrixA, int num){

  int i, j, k;
  double** result = allocate_matrix(num, num);

  for (i=num-1; i >= 0; i--){
    result[num-1][i] = matrixA[num-1][i]/matrixB[num-1][num-1];

    for (j=num-2; j>=0; j--){
      result[j][i] = matrixA[j][i];

      for (k=j+1; k < num; k++){
	result[j][i] = result[j][i] - result[k][i]*matrixB[j][k];
      }

      result[j][i] = result[j][i]/matrixB[j][j];
    }
  }
  
  return result;

}
