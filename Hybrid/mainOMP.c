#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "test_functions.h"
#include <omp.h>
#include <time.h>

#define print 0
void predictor_corrector(double**, double**, int, int, double, double);
int num_thrds;


int main(int argc, char *argv[])
{
    double **arrayA;
    int num_space, num_time, intern_grid_points;
    int  i, j;
    int row;
    int column;
    int  divisor;
    double *x;
    double *t;
    double dt;
    double dx;
    double x0 = 0.0;
    double x_end = 1.0;
    double t0 = 0.0;
    double diff_coeff = 1.0;
    double beta, alpha;
    double **solution;
 
    if (argc != 3){
      printf("Usage: <>Incomplete number of arguments, Exiting now\n");
      exit(0);
    }
    
    double t_end = atof(argv[1]);
    num_thrds = atoi(argv[2]);

    dx = 0.01;
    dt = 0.2;
    beta = 1.6;
    alpha = 0.4;
    
    num_space = (int) ((x_end - x0)/dx + 1);
    num_time = (int)  ((t_end - t0)/dt + 1);


    intern_grid_points = num_space - 2;

    arrayA = allocate_matrix(intern_grid_points, intern_grid_points);
    solution = allocate_matrix(num_time, num_space);
   
    x = allocate_vector(num_space);
    t = allocate_vector(num_time);

    x[0] = x0;
    t[0] = t0;

    for (i=1; i < num_space; i++){
      x[i] = x[i-1] + dx;
    }
    
    for (i=1; i < num_time; i++){
      t[i] = t[i-1] + dt;
    }

    arrayA = matrix_transfer(intern_grid_points, dx, beta, diff_coeff);

  
    for (i=0; i< num_space; i++)
      solution[0][i] = pow(x[i], 2) * pow((1 - x[i]),2);
    
    double t1 = timestamp();
    predictor_corrector(solution, arrayA, num_space,  num_time, alpha, dt);
    double t2 = timestamp() - t1;
    printf("%.4f\n", t2);

    
    if (! print){
      for (j= 0; j < num_space; j++){
	printf("%.3e  ",solution[num_time - 1][j]);
      }
      printf("\n");
    }
    free(x);
    free(t);
    free_matrix(solution, num_time);
    free_matrix(arrayA, intern_grid_points);
    return 0;
}

void predictor_corrector(double **Uvec, double **matrixA, int numS, int numT, double alpha, double dt){
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
  
  int i,j;
  int intern_grid_points = numS - 2;
  
  arrayB = allocate_matrix(intern_grid_points, intern_grid_points);
  identity = allocate_matrix(intern_grid_points, intern_grid_points);
  temp_mat1 = allocate_matrix(intern_grid_points, intern_grid_points);
  temp_mat2 = allocate_matrix(intern_grid_points, intern_grid_points);  

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
  arrayB = matrix_inverse(arrayB, intern_grid_points, num_thrds);

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
    //    printf("NUm-Thrds: %d\n", num_thrds);
    //    print_matrix(matrixA, intern_grid_points, intern_grid_points);


    history_term =  historyOMP(Uvec, matrixA, intern_grid_points, alpha, 0, i, i, num_thrds);

    //    print_vector(history_term, intern_grid_points);
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
