#include "test_functions.h"

void free_matrix(double  ** matrix, int num){
  int i;
  for (i=0; i < num; i++){;
    free(matrix[i]);
  }
  free(matrix);
}
