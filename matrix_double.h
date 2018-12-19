//
// Created by root on 07/12/18.
//

#ifndef ARITHMETIC_MATRIX_FACTORIZATION_H
#define ARITHMETIC_MATRIX_FACTORIZATION_H

#include <stdlib.h>
#include <stdio.h>
#include "Util.h"
#include <math.h>

#define MATRIX(matrix ,i,j) matrix.values[i * matrix.nb_col + j]
#define MATRIX_P(matrix,i,j) matrix->values[i * matrix->nb_col + j]

typedef struct matrix_double {
    unsigned int nb_line;
    unsigned int nb_col;
    double *values;
} matrix_double;


void init_matrix_double(matrix_double *matrixDouble , unsigned int nb_line , unsigned int nb_col);

void setCoeff_matrix_double(matrix_double *matrixDouble , unsigned int line , unsigned int col , double newvalue);

void setCoeff_matrix_double_array(matrix_double *matrixDouble , double *coeffs , unsigned int sizeArray);

void setAllCoeff_matrix_double(matrix_double *matrixDouble , double value);

void destroy_matrix_double(matrix_double matrixDouble);

void print_matrix_double(matrix_double matrixDouble , char *name);

void change_nb_line_matrix_double(matrix_double *matrixDouble , unsigned int new_nb_line);

void change_nb_col_matrix_double(matrix_double *matrixDouble , unsigned int new_nb_col);

void change_dim_matrix_double(matrix_double *matrixDouble , unsigned int new_nb_line , unsigned int new_nb_col);

void add_matrix_double(matrix_double *res , matrix_double A , matrix_double B);

void sub_matrix_double(matrix_double *res , matrix_double A , matrix_double B);

void scalar_mul_matrix_double(matrix_double *res , matrix_double A , double lambda);

void mul_matrix_double(matrix_double *res , matrix_double A , matrix_double B);

void scalar_div_matrix_double(matrix_double *res , matrix_double A , double lambda);

void dot_product(double *dot , matrix_double u , matrix_double v);

void copy_matrix_double(matrix_double *DEST , matrix_double SRC);

void identity_matrix_double(matrix_double *id , unsigned int size);

void matrix_line_permutation(matrix_double *P , unsigned int size , unsigned int line1 , unsigned int line2);

void matrix_col_permutation(matrix_double *Q , unsigned int size , unsigned int col1 , unsigned int col2);

void LU_decomposition_matrix_double(matrix_double *L , matrix_double *U , matrix_double A);

void swap_ligne_matrix_double(matrix_double *A , unsigned int line1 , unsigned int line2);

void swap_col_matrix_double(matrix_double *A , unsigned int col1 , unsigned int col2);

double index_max_submatrix_double(unsigned int *max_index_line , unsigned int *max_index_col , matrix_double A , unsigned int submatrix_index_line , unsigned int submatrix_index_col);

void PLUQ_decomposition(matrix_double *P , matrix_double *L , matrix_double *U , matrix_double *Q , matrix_double A);

void transpose_matrix_double(matrix_double *transpose , matrix_double A);

void setColumn_matrix_double(matrix_double *A , matrix_double vector , unsigned int indexColumn);

void getColum_matrix_double(matrix_double *column , matrix_double A , unsigned int indexColumn);

void QR_Givens(matrix_double *Q , matrix_double *R , matrix_double A);

void QR_Gram_Schimdt(matrix_double *Q , matrix_double *R , matrix_double A);

double norm_vector_double(matrix_double matrixDouble);




#endif //ARITHMETIC_MATRIX_FACTORIZATION_H
