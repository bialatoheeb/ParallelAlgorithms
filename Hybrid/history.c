#include "test_functions.h"

double *history(double **Uvec, double **matrixA, int num, double alpha, int start, int end, int col, int num_thrds){
  double *history_term;
  double coeff;
  int i, j, k;
  double temp;
  int num1 = end+1;
  history_term = allocate_vector(num);
  zero_vector(history_term, num, num_thrds);
  if (end == 0)
    return history_term;
  


#pragma omp parallel private(i, j, k, coeff,  temp) num_threads(num_thrds)  
  {
    double sum;
    double *temp_vector;
    double *nonlinear_func;
    double * temp1;
    double *temp2;
    double * temp3;
    temp_vector = allocate_vector(num);
    nonlinear_func = allocate_vector(num);
    temp1 = allocate_vector(num);
    temp2 = allocate_vector(num);
    temp3 = allocate_vector(num);
    zero_vector(temp_vector, num, 1);
    zero_vector(nonlinear_func, num, 1);
    zero_vector(temp1, num, 1);
    zero_vector(temp2, num, 1);
    zero_vector(temp3, num, 1);

    #pragma omp for schedule(dynamic) 
    for (i=start; i <= end; i++){
      for (j = 0; j < num; j++){
	temp_vector[j] = Uvec[i][j+1];
      }
     
      for (k = 0; k < num; k++){
	sum = 0.0;
	for (j = 0; j < num; j++){
	  sum += ( matrixA[k][j] * temp_vector[j] );
	}
	temp1[k] = sum;
      }
      for (k=0; k <num; k++){
	nonlinear_func[k] = pow(temp_vector[k],2);
      }
      
      
      for (k=0; k<num; k++){
	temp2[k] = temp1[k] + nonlinear_func[k];
      }
      
       
      if(i == 0)
	coeff = - (col - alpha)*pow(col+1, alpha) + pow(col, alpha)*(2*col - alpha - 1) - pow(col-1, alpha+1);
      else if (i == end && end == col)
	coeff = pow(2, alpha+1) - alpha - 3;
      else
	coeff = pow(col-i+2, alpha+1) - 3*pow(col-i+1, alpha+1) + 3*pow(col-i, alpha+1) - pow(col-i-1, alpha + 1);

      
      for (k=0; k< num; k++){
	temp1[k] = coeff* temp2[k];
      }

      //Sum locally to a vector
      for (j=0; j<num; j++)
	temp3[j] += temp1[j];
    }

    // Update the history term with each local history using a critical section
    #pragma omp critical 
     {
      for (k=0; k < num ; k++)
	history_term[k] += temp3[k];
     }
    
    free(temp1);
    free(temp2);
    free(temp3);
    free(temp_vector);
    free(nonlinear_func);
  }
  
  return history_term; 
}
