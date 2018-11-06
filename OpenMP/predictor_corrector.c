#include "test_functions.h"

void predictor_corrector(double **Uvec, double **matrixA, int numS, int numT, double alpha, double dt, int num_thrds){
  double **identity;
  double **arrayB;
  double *history_term;
  double * temp_vector;
  double *nonlinear_func1;
  double *nonlinear_func2;
  double * temp_vec1;
  double * temp_vec2;
  double * temp_vec3;
  double * temp_vec4;
  double ** temp_mat1;
  double ** temp_mat2;
  double ** lower;
  double ** upper;
  double ** pivot;
  int i,j;
  int intern_grid_points = numS - 2;
  
  arrayB = allocate_matrix(intern_grid_points, intern_grid_points);
  identity = allocate_matrix(intern_grid_points, intern_grid_points);
  temp_mat1 = allocate_matrix(intern_grid_points, intern_grid_points);
  temp_mat2 = allocate_matrix(intern_grid_points, intern_grid_points);
  
  lower = allocate_matrix(intern_grid_points, intern_grid_points);
  upper = allocate_matrix(intern_grid_points, intern_grid_points);
  pivot = allocate_matrix(intern_grid_points, intern_grid_points);
  
  history_term = allocate_vector(intern_grid_points);
  nonlinear_func1 = allocate_vector(intern_grid_points);
  nonlinear_func2 = allocate_vector(intern_grid_points);
  
  temp_vector = allocate_vector(intern_grid_points);
  temp_vec1 = allocate_vector(intern_grid_points);
  temp_vec2 = allocate_vector(intern_grid_points);
  temp_vec3 = allocate_vector(intern_grid_points);
  temp_vec4 = allocate_vector(intern_grid_points);
  
  zero_vector(history_term, intern_grid_points, num_thrds);
  identity_matrix(identity, intern_grid_points, num_thrds);
  zero_vector(temp_vector, intern_grid_points, num_thrds);
  zero_vector(nonlinear_func1, intern_grid_points, num_thrds);
  zero_vector(nonlinear_func2, intern_grid_points, num_thrds);
  zero_vector(temp_vec1, intern_grid_points, num_thrds);
  zero_vector(temp_vec2, intern_grid_points, num_thrds);
  zero_vector(temp_vec3, intern_grid_points, num_thrds);
  zero_vector(temp_vec4, intern_grid_points, num_thrds);
  

  temp_mat1 = product_matrix_scalar(identity, intern_grid_points, intern_grid_points, tgamma(alpha + 2), num_thrds);

  temp_mat2 = product_matrix_scalar(matrixA, intern_grid_points, intern_grid_points, pow(dt, alpha), num_thrds);
 

  arrayB =  add_matrix(temp_mat1, temp_mat2,intern_grid_points, intern_grid_points, num_thrds);
  arrayB = matrix_inverse(arrayB, intern_grid_points);

  temp_mat2 = product_matrix_scalar(temp_mat2, intern_grid_points, intern_grid_points, alpha, num_thrds);
  temp_mat2 = sub_matrix(temp_mat1, temp_mat2, intern_grid_points, intern_grid_points, num_thrds);

  free_matrix(temp_mat1, intern_grid_points);

  copy_matrix_vector(Uvec, temp_vector, intern_grid_points, 0, num_thrds);
  matrixA = product_matrix_scalar(matrixA, intern_grid_points, intern_grid_points, -1.0, num_thrds);

  for (i=1; i < numT; i++){    
    nonlinear_func1 = nonlinear_function(temp_vector, intern_grid_points, num_thrds);
  
  
    temp_vec1 = product_matrix_vector(temp_mat2, temp_vector, intern_grid_points, num_thrds);
    temp_vec2 = product_vector_scalar(history_term, intern_grid_points, pow(dt, alpha), num_thrds);
    
    // Predict
    temp_vec3 = product_vector_scalar(nonlinear_func1, intern_grid_points, (alpha + 1)* pow(dt, alpha), num_thrds);
    temp_vec4 = add_vector(temp_vec2, temp_vec3, intern_grid_points, num_thrds);
    temp_vec4 =  add_vector(temp_vec1, temp_vec4, intern_grid_points, num_thrds);

    temp_vec4 = product_matrix_vector(arrayB, temp_vec4, intern_grid_points, num_thrds);
   
    // Evaluate the predictor
    nonlinear_func2 = nonlinear_function(temp_vec4, intern_grid_points, num_thrds);

    //Corrector
    temp_vec3 = product_vector_scalar(nonlinear_func1, intern_grid_points, alpha*pow(dt, alpha), num_thrds);
    temp_vec4 = product_vector_scalar(nonlinear_func2, intern_grid_points, pow(dt, alpha), num_thrds);
    temp_vec4 = add_vector(add_vector(temp_vec3, temp_vec4, intern_grid_points, num_thrds), temp_vec2, intern_grid_points, num_thrds);
    temp_vec4 = add_vector(temp_vec1, temp_vec4, intern_grid_points, num_thrds);

    temp_vector = product_matrix_vector(arrayB, temp_vec4, intern_grid_points, num_thrds);
 
    copy_vector_matrix(Uvec, temp_vector, intern_grid_points, i, num_thrds);


    history_term =  historyOMP(Uvec, matrixA, intern_grid_points, alpha, 0, i, i, num_thrds);
}
  
  free_matrix(identity, intern_grid_points);
  free_matrix(arrayB, intern_grid_points);
  free_matrix(temp_mat2, intern_grid_points);
  
  free(nonlinear_func1);
  free(nonlinear_func2);
  free(history_term);
  free(temp_vector);
  free(temp_vec1);
  free(temp_vec2);
  free(temp_vec3);
}
