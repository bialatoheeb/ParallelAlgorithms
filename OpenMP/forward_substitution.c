#include <stdlib.h>
#include <stdio.h>
#include "test_functions.h"
#include <omp.h>

double ** forward_substitution(double ** matrixB, double ** matrixA, int num){

  int i, j, k;
  double** result = allocate_matrix(num, num);

  for (i=0; i<num; i++){
    result[0][i] = matrixA[0][i]/matrixB[0][0];

    for (j=1; j<num; j++){
      result[j][i] = matrixA[j][i];

      for (k=0; k < j; k++){
	result[j][i] = result[j][i] - result[k][i]*matrixB[j][k];
      }

      result[j][i] = result[j][i]/matrixB[j][j];
    }


  }
    return result;

}
