//
// Created by root on 29/10/18.
//

#ifndef ARITHMETIC_POL_BIGFLOAT_ARITHMETIC_H
#define ARITHMETIC_POL_BIGFLOAT_ARITHMETIC_H


#include <gmp.h>
#include "Util.h"
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>

typedef struct pol_bigfloat {
    unsigned int degree;
    mpf_t *coeffs;
} pol_bigfloat;


void init_pol_bigfloat(pol_bigfloat *polynomial, unsigned int degree);

void set_all_coeffs_pol_bigfloat(pol_bigfloat polynomial, mpf_t *coeffs, size_t size_coeffs_array);

void set_all_coeffs_random_pol_bigfloat(pol_bigfloat polynomial, unsigned int max);

void set_coeff_pol_bigfloat(pol_bigfloat polynomial, unsigned int degree_coeff, mpf_t newValue);

void set_all_coeffs_to_pol_bigfloat(pol_bigfloat polynomial, mpf_t value);

void destroy_pol_bigfloat(pol_bigfloat polynomial);

void change_degre_pol_bigfloat(pol_bigfloat *polynomial, unsigned int new_degree);

void copy_pol_bigfloat(pol_bigfloat *res, pol_bigfloat polynomial);

void print_pol_bigfloat(pol_bigfloat polynomial);

void add_pol_bigfloat(pol_bigfloat *res, pol_bigfloat A, pol_bigfloat B);

void sub_pol_bigfloat(pol_bigfloat *res, pol_bigfloat A, pol_bigfloat B);

void mult_pol_bigfloat(pol_bigfloat *res, pol_bigfloat A, pol_bigfloat B);

void karatsuba_pol_bigfloat(pol_bigfloat *res, pol_bigfloat A, pol_bigfloat B);

#endif //ARITHMETIC_POL_BIGFLOAT_ARITHMETIC_H
