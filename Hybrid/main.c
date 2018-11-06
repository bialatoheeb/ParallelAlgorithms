#include "test_functions.h"
#include "mpi.h"

#define MAX_RANKS 64
#define ROOT 0
#define PRINT 1 // Change to 1 if you wish to print the solution at final time T

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
    int num_ranks;
    int my_rank;
    
    if (argc != 3){
      printf("Usage: <>Incomplete number of arguments:  ");
      printf("mpiexec -n <num_ranks> ./main <Final Time> <# of threads>.\nExiting now\n");
      exit(0);
    }
     
    double t_end = atof(argv[1]);
    int num_thrds = atoi(argv[2]);
    dx = 0.05;
    dt = 0.1;
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
    
    for (i=1; i< num_time; i++){
      t[i] = t[i-1] + dt;
    }

    arrayA = matrix_transfer(intern_grid_points, dx, beta, diff_coeff);
   
    for (i=0; i< num_space; i++)
      solution[0][i] = pow(x[i], 2) * pow((1 - x[i]),2);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    double t1 = timestamp();
    predictor_corrector(solution, arrayA, num_space,  num_time, alpha, dt, num_thrds);
    double t2 = timestamp() - t1;
    if (my_rank == ROOT){
      printf("%.3f\n", t2);
      
      if (PRINT){
	for (j= 0; j < num_space; j++){
	  printf("%.3e  ", solution[num_time - 1][j]);
	}
	
	printf("\n");
      }
    }
    
    free(x);
    free(t);
    free_matrix(solution, num_time);
    free_matrix(arrayA, intern_grid_points);    
    MPI_Finalize();
    return 0;
}

