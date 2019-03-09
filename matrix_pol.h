//
// Created by micky on 27/02/19.
//

#ifndef ARITHMETIC_MATRIX_pol_POL_H
#define ARITHMETIC_MATRIX_pol_POL_H

#include "pol.h"


#define MATRIX(matrix ,i,j) matrix.values[i * matrix.nb_col + j]
#define MATRIX_P(matrix,i,j) matrix->values[i * matrix->nb_col + j]


typedef struct matrix_pol {
    unsigned int nb_line;
    unsigned int nb_col;
    pol *values;
} matrix_pol;

void init_matrix_pol(matrix_pol *matrix_pol , unsigned int nb_line , unsigned int nb_col);

void setCoeff_matrix_pol(matrix_pol *matrix_pol , unsigned int line , unsigned int col , pol newvalue);

void set_allCoeff_matrix_pol(matrix_pol *matrix_pol , pol newvalue);

void destroy_matrix_pol(matrix_pol matrix_pol);

void print_matrix_pol(matrix_pol matrix_pol , char *name);

void change_nb_line_matrix_pol(matrix_pol *matrix_pol , unsigned int new_nb_line);

void change_nb_col_matrix_pol(matrix_pol *matrix_pol , unsigned int new_nb_col);

void change_dim_matrix_pol(matrix_pol *matrix_pol , unsigned int new_nb_line , unsigned int new_nb_col);

void add_matrix_pol(matrix_pol *res , matrix_pol A , matrix_pol B);

void mul_matrix_pol(matrix_pol *res , matrix_pol A , matrix_pol B);

void copy_matrix_pol(matrix_pol *DEST , matrix_pol SRC);

void identity_matrix_pol(matrix_pol *id, unsigned int size);


#endif //ARITHMETIC_MATRIX_pol_POL_H
