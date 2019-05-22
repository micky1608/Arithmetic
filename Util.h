//
// Created by root on 13/10/18.
//

#ifndef ARITHMETIC_UTIL_H
#define ARITHMETIC_UTIL_H

#include <math.h>
#include <gmp.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define MATRIX(matrix ,i,j) matrix.values[i * matrix.nb_col + j]
#define MATRIX_P(matrix,i,j) matrix->values[i * matrix->nb_col + j]

typedef enum bool { TRUE, FALSE } bool;

int max(int a , int b);

double log_base_2(double x);

void mpq_mul_si(mpq_t res , mpq_t x , int lambda_num , unsigned int lambda_den);

void mpq_pow(mpq_t res , mpq_t x , int p);

void swap(double *a , double *b);

void swap_long(long *a , long *b);

typedef struct pol {
    unsigned int degree;
    long *coeffs;
} pol;

typedef struct matrix {
    unsigned int nb_line;
    unsigned int nb_col;
    long *values;
} matrix;


#endif //ARITHMETIC_UTIL_H
