//
// Created by root on 12/10/18.
//

#ifndef ARITHMETIC_POLYNOMIAL_ARITHMETIC_H
#define ARITHMETIC_POLYNOMIAL_ARITHMETIC_H

#include <gmp.h>
#include "Util.h"
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>

typedef struct pol_bigint {
    unsigned int degree;
    mpz_t *coeffs;
} pol_bigint;


void init_pol_bigint(pol_bigint *polynomial, unsigned int degree);

void set_all_coeffs_pol_bigint(pol_bigint polynomial, mpz_t *coeffs, size_t size_coeffs_array);

void set_all_coeffs_random_pol_bigint(pol_bigint polynomial, unsigned int max);

void set_coeff_pol_bigint(pol_bigint polynomial, unsigned int degree_coeff, mpz_t newValue);

void set_all_coeffs_to_pol_bigint(pol_bigint polynomial, mpz_t value);

void destroy_pol_bigint(pol_bigint polynomial);

void change_degre_pol_bigint(pol_bigint *polynomial, unsigned int new_degree);

void copy_pol_bigint(pol_bigint *res, pol_bigint polynomial);

void print_pol_bigint(pol_bigint polynomial);

void add_pol_bigint(pol_bigint *res, pol_bigint A, pol_bigint B);

void sub_pol_bigint(pol_bigint *res, pol_bigint A, pol_bigint B);

void mult_pol_bigint(pol_bigint *res, pol_bigint A, pol_bigint B);

void karatsuba_pol_bigint(pol_bigint *res, pol_bigint A, pol_bigint B);


#endif //ARITHMETIC_POLYNOMIAL_ARITHMETIC_H
