CC= gcc -fopenmp
HEADER = test_functions.h
FLAG = #-O2
# -I includes all depndent files

PROGS = add_matrix.o\
	add_vector.o\
	timestamp.o\
	permutation_matrix.o\
	LU_decomposition.o\
	forward_substitution.o\
	backward_substitution.o\
	matrix_inverse.o\
	allocate_matrix.o\
	allocate_vector.o\
	copy_matrix.o\
	copy_vector_matrix.o\
	copy_matrix_vector.o\
	free_matrix.o\
	historyOMP.o\
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



mainOMP: mainOMP.o $(PROGS)
	$(CC) -lm -o mainOMP mainOMP.o $(PROGS) -I. $(FLAG)

mainOMP.o: mainOMP.c
	$(CC) -c -o mainOMP.o mainOMP.c -I. $(FLAG)

%.o: %.c $(HEADER)
	$(CC) -c -o $@ $< -I. $(FLAG)


.PHONY:


clean:
	rm -f *.o *~  mainOMP makefile~
