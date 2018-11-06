#include "test_functions.h"

double **matrix_transfer( int num, double dx, double alpha, double diff_coeff){
  int i,j, m, n;
  double product;
  double **temp1 = allocate_matrix(num, num);
  double **temp2 = allocate_matrix(num, num);
  double **temp3 = allocate_matrix(num, num);
  double **inverse;
  double **array;


  zero_matrix(temp2, num, num);

  for(i=1; i <= num; i++){
    for (j=1; j <= num; j++){
      temp2[i-1][j-1] = sin(PI*i*j/(num+1));
    }
    product = 2.0*sin(i*PI/(2.0*(num+1)));
    temp1[i-1][i-1] = diff_coeff/pow(dx, alpha)*pow(product, alpha);
  }
  inverse = copy_matrix(temp2, num, num);

  inverse = matrix_inverse(inverse, num);

  temp3 = product_matrix_matrix(temp1, temp2, num);
  array = product_matrix_matrix(inverse, temp3, num);

  free_matrix(temp1, num);
  free_matrix(temp2, num);
  free_matrix(temp3, num);
  free_matrix(inverse, num);

  return array;
}
