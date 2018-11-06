CC= gcc
CMPI=mpicc
CFLAGS=
HEADER = test_functions.h


PROGS = add_matrix.o\
	add_vector.o\
	timestamp.o\
	backward_substitution.o\
	forward_substitution.o\
	LU_decomposition.o\
	matrix_inverse.o\
	permutation_matrix.o\
	allocate_matrix.o\
	allocate_vector.o\
	copy_matrix.o\
	copy_vector_matrix.o\
	copy_matrix_vector.o\
	free_matrix.o\
	history.o\
	identity_matrix.o\
	matrix_transfer.o\
	nonlinear_function.o\
	predictor_corrector.o\
	print_matrix.o\
	print_vector.o\
	product_matrix_matrix.o\
	product_matrix_scalar.o\
	product_matrix_vector.o\
	product_vector_scalar.o\
	sub_matrix.o\
	zero_matrix.o\
	zero_vector.o




main: main.o $(PROGS)
	$(CMPI) -lm -o main main.o $(PROGS) -I.

main.o: main.c 
	$(CMPI) -c main.c -I.

predictor_corrector.o: predictor_corrector.c  
	$(CMPI) -lm -c predictor_corrector.c -I.

%.o: %.c $(HEADER)
	$(CC) -c -o $@ $< -I.


.PHONY:


clean:
	rm -f *.o *~  main 
