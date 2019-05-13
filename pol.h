#ifndef ARITHMETIC_POLYNOMIAL_ARITHMETIC_H
#define ARITHMETIC_POLYNOMIAL_ARITHMETIC_H

#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include "arithmetic.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef struct pol {
    unsigned int degree;
    long *coeffs;
} pol;


void init_pol(pol *polynomial, unsigned int degree);

void set_all_coeffs_pol(pol polynomial, long *coeffs, size_t size_coeffs_array);

void set_all_coeffs_random_pol(pol polynomial, unsigned int max);

void set_coeff_pol(pol polynomial, unsigned int degree_coeff, long newValue);

void set_all_coeffs_to_pol(pol polynomial, long value);

void destroy_pol(pol polynomial);

void change_degre_pol(pol *polynomial, unsigned int new_degree);

void copy_pol(pol *res, pol polynomial);

void print_pol(pol polynomial , char *name);

void print_pol_center(pol polynomial , char *name , unsigned int size);

void add_pol(pol *res, pol A, pol B);

void sub_pol(pol *res, pol A, pol B);

void mult_pol(pol *res, pol A, pol B);

void scalar_mult_pol(pol *res , pol A , long lambda);

int is_zero_pol(pol polynomial);

void euclide_div_pol(pol *Q , pol *R , pol A , pol B);

void horner_eval(long *res , pol f , long x);

void horner_eval_multi(long *res , pol f , long *x , unsigned int nb_x);

void derivate_pol(pol *derivate , pol A);

/* 'ff' stands for finite field. I generally work in the fields Z/pZ, namely {0,1,2 ... p-1} with p a prime number. */

void reduce_pol_ff(pol *A , long P); 

void euclide_div_pol_ff(pol *Q , pol *R , pol A , pol B , long P);

void gcd_pol_ff(pol *gcd , pol A , pol B , long P); 


#endif //ARITHMETIC_POLYNOMIAL_ARITHMETIC_H