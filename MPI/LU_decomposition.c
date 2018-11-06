#include "test_functions.h"

void LU_decomposition(double **matrix, double **L, double **U, double** pivot, int num){
  int i,j, k;
  double **array = allocate_matrix(num, num);
   
  zero_matrix(L, num, num);
  zero_matrix(U, num, num);
  permutation_matrix(matrix, pivot, num);

  array =  product_matrix_matrix(pivot, matrix, num);

  for(i=0; i < num; i++){
    L[i][i] = 1.0;
  }
 
  for (i=0; i < num; i++){
    for (j=0; j < num; j++){
      double temp;
      if (j <=i){
	temp =0;
	for (k=0; k < j; k++)
	  temp += (L[j][k]*U[k][i]);
	U[j][i] = array[j][i] - temp;
      }
      if (j >= i){
	temp = 0;
	for (k=0; k < i; k++)
	  temp += (L[j][k]*U[k][i]);
	L[j][i] = (array[j][i] - temp)/U[i][i];
      }
    }
  }

  free_matrix(array, num);
}
