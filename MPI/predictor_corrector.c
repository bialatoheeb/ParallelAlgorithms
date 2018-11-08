#include "test_functions.h"
#include "mpi.h"

#define ROOT 0
void  predictor_corrector(double **Uvec, double **matrixA, int numS, int numT, double alpha, double dt){
  double **identity;
  double **arrayB;
  double *history_term1;
  double *history_term2;
  
  double * temp_vector;
  double * temp_array;
  double *nonlinear_func1;
  double *nonlinear_func2;
  double * temp_vec1;
  double * temp_vec2;
  double * temp_vec3;
  double *bsend_array;
  double ** temp_mat2;
  double dt_alpha;
  int my_rank, num_ranks;
  int i,j, k;
  int num_blocks; 
  int intern_grid_points = numS - 2;
  int buffer, val, tag;
  MPI_Status status;
  MPI_Datatype send_type, broadcast_type;

  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  arrayB = allocate_matrix(intern_grid_points, intern_grid_points);
  identity = allocate_matrix(intern_grid_points, intern_grid_points);
  temp_mat2 = allocate_matrix(intern_grid_points, intern_grid_points);
  bsend_array = allocate_vector(num_ranks* numS);
  temp_array = allocate_vector((my_rank+1)*numS);
  history_term1 = allocate_vector(intern_grid_points);
  history_term2 = allocate_vector(intern_grid_points);
  

  nonlinear_func1 = allocate_vector(intern_grid_points);
  nonlinear_func2 = allocate_vector(intern_grid_points);
  
  temp_vector = allocate_vector(intern_grid_points);
  temp_vec1 = allocate_vector(intern_grid_points);
  temp_vec2 = allocate_vector(intern_grid_points);
  temp_vec3 = allocate_vector(intern_grid_points);
    
  zero_vector(history_term1, intern_grid_points);
  zero_vector(history_term2, intern_grid_points);
  zero_vector(temp_vector, intern_grid_points);
  zero_vector(nonlinear_func1, intern_grid_points);
  zero_vector(nonlinear_func2, intern_grid_points);
  zero_vector(temp_vec1, intern_grid_points);
  zero_vector(temp_vec2, intern_grid_points);
  zero_vector(temp_vec3, intern_grid_points);
  zero_vector(bsend_array, num_ranks* numS);
  identity_matrix(identity, intern_grid_points);

  MPI_Type_contiguous(numS*(my_rank + 1), MPI_DOUBLE, &send_type);
  MPI_Type_contiguous(numS*num_ranks, MPI_DOUBLE, &broadcast_type);
  MPI_Type_commit(&send_type);
  MPI_Type_commit(&broadcast_type);
  
  dt_alpha = pow(dt, alpha);
  num_blocks = (int)((numT - 1)/num_ranks);

  // Precompute matrix P = arrayB and Q = matrixA  as described in the paper.
  temp_mat2 = product_matrix_scalar(identity, intern_grid_points, intern_grid_points, tgamma(alpha + 2));

  arrayB = add_matrix(
		      temp_mat2,
		      product_matrix_scalar(matrixA, intern_grid_points, intern_grid_points, dt_alpha),
		      intern_grid_points, intern_grid_points
		      );
  arrayB = matrix_inverse(arrayB, intern_grid_points);

  temp_mat2 = sub_matrix(
			 temp_mat2,
			 product_matrix_scalar(matrixA, intern_grid_points, intern_grid_points, alpha*dt_alpha),
			 intern_grid_points, intern_grid_points
			 );
  matrixA = product_matrix_scalar(matrixA, intern_grid_points, intern_grid_points, -1.0);
 
  val = 0;  
  for (i = 0; i < num_blocks; i++){

    // At  the beginning of each iteration each rank copies the most recent solution
    // to its array of solution
    for (k=1; k <= num_ranks; k++){
      for (j=0; j < numS; j++){
	Uvec[val+k][j] = bsend_array[(k-1)*numS + j];
      }
    }

    // The next two line is the parallel part of this algroithm
    // Each rank computes I3 given in the paper
    val = i*num_ranks;
    history_term1 = history(Uvec, matrixA, intern_grid_points, alpha, 0, val, val + my_rank);

    tag = my_rank - 1;
    zero_vector(temp_array, (my_rank+1) *numS);

    // Recv result from rank preceding the present rank except for the root
    if (my_rank != ROOT)
      MPI_Recv(temp_array, 1, send_type, tag, tag, MPI_COMM_WORLD, &status);

    // After recv, copy to array of solution for use with the history term
    for (k=1; k <= my_rank; k++){
      for (j=0; j < numS; j++){
	Uvec[val + k][j] = temp_array[(k-1)*numS + j];
      }
    }

   
    // All ranks except ROOT computes the other "half" of thier history term 
    if (my_rank != ROOT){
      if (i != 0)
	history_term2 = history(Uvec, matrixA, intern_grid_points, alpha, val + 1, val + my_rank, val + my_rank);
      else
	history_term2 = history(Uvec, matrixA, intern_grid_points, alpha, 0, val + my_rank, val + my_rank);
    }

    // Add up both history terms
    history_term1 = add_vector(history_term1, history_term2, intern_grid_points);

    // Get the vector at  time i-1 if I'm computing at time i;
    copy_matrix_vector(Uvec, temp_vector, intern_grid_points, val + my_rank);

    // Compute the nonlinear function
    nonlinear_func1 = nonlinear_function(temp_vector, intern_grid_points);
    
    temp_vec1 = product_matrix_vector(temp_mat2, temp_vector, intern_grid_points);
    temp_vec2 = product_vector_scalar(history_term1, intern_grid_points, dt_alpha);
    temp_vec2 = add_vector(temp_vec1, temp_vec2, intern_grid_points);
      

    // Predict
    temp_vec3 = add_vector(
			   temp_vec2,
			   product_vector_scalar(nonlinear_func1, intern_grid_points, (alpha + 1)* dt_alpha),
			   intern_grid_points
			   );
    // Predictor
    temp_vec3 = product_matrix_vector(arrayB, temp_vec3, intern_grid_points);

    // Evaluate the nonlinear function using the predictor
    nonlinear_func2 = nonlinear_function(temp_vec3, intern_grid_points);

    //Correct
    temp_vec3 = add_vector(
			   add_vector(
				      product_vector_scalar(nonlinear_func1, intern_grid_points, alpha*dt_alpha),
				      product_vector_scalar(nonlinear_func2, intern_grid_points, dt_alpha),
				      intern_grid_points
				      ),
			   temp_vec2,
			   intern_grid_points
			   );
    // Corrector
    temp_vector = product_matrix_vector(arrayB, temp_vec3, intern_grid_points);      

    // Copy new solution to temp_array
    for (k = 0; k < intern_grid_points; k++)
      temp_array[(my_rank*numS) + k + 1] =  temp_vector[k];

    // Send (to my my_rank + 1) only if I'm not the last rank, 
    // otherwise copy to the bsend_array for broadcast to all ranks
    tag = my_rank;
    if (my_rank != (num_ranks - 1)){
      MPI_Send( temp_array, 1, send_type, my_rank+1, tag, MPI_COMM_WORLD);
    }else{
      for (k=0; k < num_ranks; k++){
	for (j=0; j< numS; j++){
	  bsend_array[k*numS + j] = temp_array[k*numS + j];
	}
      }
    }

    MPI_Bcast(bsend_array, 1, broadcast_type, num_ranks - 1, MPI_COMM_WORLD );
  }

  // Copy the last broadcast solution to the array of solution
  for (k=1; k <= num_ranks; k++){
    for (j=0; j < numS; j++){
      Uvec[val+k][j] = bsend_array[(k-1)*numS + j];
    }
  }
  
  // Check for remainder of the t-stencil and perform one iteration
  // of the for-loop above if remainder is not 0
  
  int rem = (numT - 1) % num_ranks;
  if (rem != 0){
    val = num_ranks * num_blocks;
    if ( my_rank < rem){
    history_term1 = history(Uvec, matrixA, intern_grid_points, alpha, 0, val, val + my_rank);

    tag = my_rank - 1;
    zero_vector(temp_array, (my_rank+1) *numS);

    if (my_rank != ROOT)
      MPI_Recv(temp_array, 1, send_type, tag, tag, MPI_COMM_WORLD, &status);


    for (k=1; k <= my_rank; k++){
      for (j=0; j < numS; j++){
	Uvec[val + k][j] = temp_array[(k-1)*numS + j];
      }
    }

    if (my_rank != ROOT){
      if (i != 0)
	history_term2 = history(Uvec, matrixA, intern_grid_points, alpha, val + 1, val + my_rank, val + my_rank);
      else
	history_term2 = history(Uvec, matrixA, intern_grid_points, alpha, 0, val + my_rank, val + my_rank);
    }
    history_term1 = add_vector(history_term1, history_term2, intern_grid_points);

    copy_matrix_vector(Uvec, temp_vector, intern_grid_points, val + my_rank);

    nonlinear_func1 = nonlinear_function(temp_vector, intern_grid_points);

    temp_vec1 = product_matrix_vector(temp_mat2, temp_vector, intern_grid_points);
    temp_vec2 = product_vector_scalar(history_term1, intern_grid_points, dt_alpha);
    temp_vec2 = add_vector(temp_vec1, temp_vec2, intern_grid_points);


    // Predict
    temp_vec3 = add_vector(
			   temp_vec2,
			   product_vector_scalar(nonlinear_func1, intern_grid_points, (alpha + 1)* dt_alpha),
			   intern_grid_points
			   );
    // Predictor
    temp_vec3 = product_matrix_vector(arrayB, temp_vec3, intern_grid_points);

    // Evaluate nonlinear function with the predictor
    nonlinear_func2 = nonlinear_function(temp_vec3, intern_grid_points);

    //Correct
    temp_vec3 = add_vector(
			   add_vector(
				      product_vector_scalar(nonlinear_func1, intern_grid_points, alpha*dt_alpha),
				      product_vector_scalar(nonlinear_func2, intern_grid_points, dt_alpha),
				      intern_grid_points
				      ),
			   temp_vec2,
			   intern_grid_points
			   );
    // Corrector
    temp_vector = product_matrix_vector(arrayB, temp_vec3, intern_grid_points);

    for (k = 0; k < intern_grid_points; k++)
      temp_array[(my_rank*numS) + k + 1] =  temp_vector[k];

    tag = my_rank;
    if (my_rank <  (rem -1)){
      MPI_Send( temp_array, 1, send_type, my_rank+1, tag, MPI_COMM_WORLD);
    }else{
      for (k=0; k < rem; k++){
	for (j=0; j< numS; j++){
	  bsend_array[k*numS + j] = temp_array[k*numS + j];
	}
      }
    }
    }
    MPI_Bcast(bsend_array, 1, broadcast_type, rem - 1, MPI_COMM_WORLD );

    for (k=1; k <= rem; k++){
      for (j=0; j < numS; j++){
	Uvec[val+k][j] = bsend_array[(k-1)*numS + j];
      }
    }
  }

  // Free all allocated vectors and arrays  
  free_matrix(identity, intern_grid_points);
  free_matrix(arrayB, intern_grid_points);
  free_matrix(temp_mat2, intern_grid_points);
  free(bsend_array);
  free(temp_array);
  free(nonlinear_func1);
  free(nonlinear_func2);
  free(history_term1);
  free(history_term2);
  free(temp_vector);
  free(temp_vec1);
  free(temp_vec2);
  free(temp_vec3);
}


