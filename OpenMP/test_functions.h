#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#define PI 3.14159265358979323846


double timestamp();
void predictor_corrector(double **, double **, int, int, double, double, int);
double * nonlinear_function(double *, int, int);                                                                                                            
double** matrix_transfer(int, double, double, double);                                                                                                
void inverse_matrix(double **, int , int);
void zero_matrix(double**,int,int, int);                                                                                                                    
void zero_vector(double * , int, int);                                                                                                                           
void free_matrix(double **, int);                                                                                                            
double **allocate_matrix(int, int);                                                                                                                         
double *allocate_vector(int);                                                                                                                               
void identity_matrix(double **, int, int);                                                                                                                       
double **product_matrix_matrix(double **, double **, int, int);                                                                                                  
double *product_matrix_vector(double **, double *, int, int);                                                                                                    
double **product_matrix_scalar(double **, int, int, double, int);                                                                                                
double *product_vector_scalar(double *, int, double, int);                                                                                                       
void copy_matrix_vector(double **, double *, int, int, int);                                                                                                     
void copy_vector_matrix(double **, double *, int, int, int);                                                                                                     
double **copy_matrix(double **, int, int, int);                                                                                                                  
double **add_matrix(double **, double **, int, int, int);                                                                                                        
double **sub_matrix(double **, double **, int, int, int);                                                                                                        
double *add_vector(double *, double *, int, int);                                                                                                                
double *sub_vector(double *, double *, int, int);                                                                                                             
double *historyOMP(double **, double **, int, double, int, int, int, int);
double *add_vector(double *, double *, int, int);
double *add_OMP(double *, double *, int, int, int);
void print_matrix(double **, int, int);
void print_vector(double *, int);
void permutation_matrix(double **, double **, int);
void LU_decomposition(double **, double **, double **, double **, int);
double ** forward_substitution(double **, double **, int);
double ** backward_substitution(double **, double **, int);
double * forward(double **, double **, double *, int, int);
double * backward(double **, double *, int, int);
double ** matrix_inverse(double **, int);
