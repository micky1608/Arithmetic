//
// Created by root on 30/10/18.
//

#ifndef ARITHMETIC_MATRIX_BIGINT_ARITHMETIC_H
#define ARITHMETIC_MATRIX_BIGINT_ARITHMETIC_H

#include <gmp.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

typedef struct matrix_bigint {
    unsigned int nb_line;
    unsigned int nb_col;
    mpz_t *values;
} matrix_bigint;

void init_matrix_bigint(matrix_bigint *matrix , unsigned int nb_line_p , unsigned nb_col_p);

void set_coeff_matrix_bigint(matrix_bigint matrix, unsigned int index_line , unsigned int index_col , mpz_t newvalue);

void set_coeff_matrix_bigint_d(matrix_bigint matrix, unsigned int index_line , unsigned int index_col , int newvalue);

void set_all_coeffs_random_matrix_bigint(matrix_bigint matrix , unsigned int max);

void destroy_matrix_bigint(matrix_bigint matrix);

void print_matrix_bigint(matrix_bigint matrix);

void change_nb_line_matrix_bigint(matrix_bigint *matrix , unsigned int new_nb_line);

void change_nb_col_matrix_bigint(matrix_bigint *matrix , unsigned int new_nb_col);

void change_dim_matrix_bigint(matrix_bigint *matrix , unsigned int new_nb_line , unsigned int new_nb_col);

#endif //ARITHMETIC_MATRIX_BIGINT_ARITHMETIC_H
