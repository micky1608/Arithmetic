//
// Created by root on 30/10/18.
//

#ifndef ARITHMETIC_MATRIX_BIGQ_ARITHMETIC_H
#define ARITHMETIC_MATRIX_BIGQ_ARITHMETIC_H

#include <gmp.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

typedef struct matrix_bigQ {
    unsigned int nb_line;
    unsigned int nb_col;
    mpq_t *values;
} matrix_bigQ;

void init_matrix_bigQ(matrix_bigQ *matrix , unsigned int nb_line_p , unsigned nb_col_p);

void set_coeff_matrix_bigQ(matrix_bigQ matrix, unsigned int index_line , unsigned int index_col , mpq_t newvalue);

void set_coeff_matrix_bigQ_d(matrix_bigQ matrix, unsigned int index_line , unsigned int index_col , int numerator , unsigned int denominator);

void set_all_coeffs_random_matrix_bigQ(matrix_bigQ matrix , unsigned int max);

void destroy_matrix_bigQ(matrix_bigQ matrix);

void print_matrix_bigQ(matrix_bigQ matrix);

void change_nb_line_matrix_bigQ(matrix_bigQ *matrix , unsigned int new_nb_line);

void change_nb_col_matrix_bigQ(matrix_bigQ *matrix , unsigned int new_nb_col);

void change_dim_matrix_bigQ(matrix_bigQ *matrix , unsigned int new_nb_line , unsigned int new_nb_col);

void add_matrix_bigQ(matrix_bigQ *res , matrix_bigQ A , matrix_bigQ B);

void sub_matrix_bigQ(matrix_bigQ *res , matrix_bigQ A , matrix_bigQ B);

void scalar_mult_matrix_bigQ(matrix_bigQ *res , matrix_bigQ A , long lambda);

void mult_matrix_bigQ(matrix_bigQ *res , matrix_bigQ A , matrix_bigQ B);

#endif //ARITHMETIC_MATRIX_BIGQ_ARITHMETIC_H
