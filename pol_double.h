#ifndef ARITHMETIC_POLYNOMIAL_ARITHMETIC_DOUBLE_H
#define ARITHMETIC_POLYNOMIAL_ARITHMETIC_DOUBLE_H

#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef struct {
    unsigned int degree;
    double *coeffs;
} pol_double;


void init_pol_double(pol_double* polynomial, unsigned int degree);

void set_all_coeffs_pol_double(pol_double polynomial, double *coeffs, size_t size_coeffs_array);

void set_all_coeffs_random_pol_double(pol_double polynomial, unsigned int max);

void set_coeff_pol_double(pol_double polynomial, unsigned int degree_coeff, double newValue);

void set_all_coeffs_to_pol_double(pol_double polynomial, double value);

void destroy_pol_double(pol_double polynomial);

void change_degre_pol_double(pol_double *polynomial, unsigned int new_degree);

void copy_pol_double(pol_double *res, pol_double polynomial);

void print_pol_double(pol_double polynomial , char *name);

void print_pol_double_center(pol_double polynomial , char *name , unsigned int size);

void add_pol_double(pol_double *res, pol_double A, pol_double B);

void sub_pol_double(pol_double *res, pol_double A, pol_double B);

void mult_pol_double(pol_double *res, pol_double A, pol_double B);

void scalar_mult_pol_double(pol_double *res , pol_double A , double lambda);

int is_zero_pol_double(pol_double polynomial);

void euclide_div_pol_double(pol_double *Q , pol_double *R , pol_double A , pol_double B);

void horner_eval_double(double *res , pol_double f , double x);

void horner_eval_multi_double(double *res , pol_double f , double *x , unsigned int nb_x);


#endif //ARITHMETIC_pol_doubleYNOMIAL_ARITHMETIC_DOUBLE_H