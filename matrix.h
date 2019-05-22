//
// Created by root on 07/12/18.
//

#ifndef ARITHMETIC_MATRIX_LONG_FACTORIZATION_H
#define ARITHMETIC_MATRIX_LONG_FACTORIZATION_H

#include <stdlib.h>
#include <stdio.h>
#include "Util.h"
#include "arithmetic.h"
#include <math.h>


void init_matrix(matrix *matrix , unsigned int nb_line , unsigned int nb_col);

void setCoeff_matrix(matrix *matrix , unsigned int line , unsigned int col , long newvalue);

void setCoeff_matrix_array(matrix *matrix , long *coeffs , unsigned int sizeArray);

void setAllCoeff_matrix(matrix *matrix , long value);

void destroy_matrix(matrix matrix);

void print_matrix(matrix matrix , char *name);

void change_nb_line_matrix(matrix *matrix , unsigned int new_nb_line);

void change_nb_col_matrix(matrix *matrix , unsigned int new_nb_col);

void change_dim_matrix(matrix *matrix , unsigned int new_nb_line , unsigned int new_nb_col);

void add_matrix(matrix *res , matrix A , matrix B);

void sub_matrix(matrix *res , matrix A , matrix B);

void scalar_mul_matrix(matrix *res , matrix A , long lambda);

void mul_matrix(matrix *res , matrix A , matrix B);

void scalar_div_matrix(matrix *res , matrix A , long lambda);

void dot_product(long *dot , matrix u , matrix v);

void copy_matrix(matrix *DEST , matrix SRC);

void identity_matrix(matrix *id , unsigned int size);

void matrix_line_permutation(matrix *P , unsigned int size , unsigned int line1 , unsigned int line2);

void matrix_col_permutation(matrix *Q , unsigned int size , unsigned int col1 , unsigned int col2);

void LU_decomposition_matrix(matrix *L , matrix *U , matrix A);

void LU_decomposition_matrix_ff(matrix *L , matrix *U , matrix A, long p);

void swap_ligne_matrix(matrix *A , unsigned int line1 , unsigned int line2);

void swap_col_matrix(matrix *A , unsigned int col1 , unsigned int col2);

long index_max_submatrix(unsigned int *max_index_line , unsigned int *max_index_col , matrix A , unsigned int submatrix_index_line , unsigned int submatrix_index_col);

void transpose_matrix(matrix *transpose , matrix A);

void setColumn_matrix(matrix *A , matrix vector , unsigned int indexColumn);

void getColum_matrix(matrix *column , matrix A , unsigned int indexColumn);

long norm_vector_long(matrix matrix);

void sylvester_matrix(matrix *syl , pol *A , pol *B);

void determinant_matrix_ff(long *det , matrix A, long p);

void reduce_matrix_ff(matrix A , long p);


#endif //ARITHMETIC_MATRIX_LONG_FACTORIZATION_H
