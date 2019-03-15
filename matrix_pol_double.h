//
// Created by micky on 27/02/19.
//

#ifndef ARITHMETIC_MATRIX_POL_DOUBLE_H
#define ARITHMETIC_MATRIX_POL_DOUBLE_H

#include "pol_double.h"


#define MATRIX(matrix ,i,j) matrix.values[i * matrix.nb_col + j]
#define MATRIX_P(matrix,i,j) matrix->values[i * matrix->nb_col + j]


typedef struct {
    unsigned int nb_line;
    unsigned int nb_col;
    pol_double *values;
} matrix_pol_double;

void init_matrix_pol_double(matrix_pol_double *matrix_pol_double , unsigned int nb_line , unsigned int nb_col);

void setCoeff_matrix_pol_double(matrix_pol_double *matrix_pol_double , unsigned int line , unsigned int col , pol_double newvalue);

void set_allCoeff_matrix_pol_double(matrix_pol_double *matrix_pol_double , pol_double newvalue);

void set_coeff_constant_matrix_pol_double(matrix_pol_double *matrix_pol_double , unsigned int line , unsigned int col , long value);

void destroy_matrix_pol_double(matrix_pol_double matrix_pol_double);

void print_matrix_pol_double(matrix_pol_double matrix_pol_double , char *name);

void change_nb_line_matrix_pol_double(matrix_pol_double *matrix_pol_double , unsigned int new_nb_line);

void change_nb_col_matrix_pol_double(matrix_pol_double *matrix_pol_double , unsigned int new_nb_col);

void change_dim_matrix_pol_double(matrix_pol_double *matrix_pol_double , unsigned int new_nb_line , unsigned int new_nb_col);

void add_matrix_pol_double(matrix_pol_double *res , matrix_pol_double A , matrix_pol_double B);

void mul_matrix_pol_double(matrix_pol_double *res , matrix_pol_double A , matrix_pol_double B);

void copy_matrix_pol_double(matrix_pol_double *DEST , matrix_pol_double SRC);

void identity_matrix_pol_double(matrix_pol_double *id, unsigned int size);


#endif // ARITHMETIC_MATRIX_POL_DOUBLE_H
