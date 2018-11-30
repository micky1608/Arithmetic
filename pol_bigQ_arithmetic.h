//
// Created by root on 30/11/18.
//

#ifndef ARITHMETIC_POL_BIGQ_ARITHMETIC_H
#define ARITHMETIC_POL_BIGQ_ARITHMETIC_H

#include <gmp.h>

typedef struct pol_bigQ {
    unsigned int degree;
    mpq_t *coeffs;
} pol_bigQ;


void init_pol_bigQ(pol_bigQ *polynomial, unsigned int degree);

void set_all_coeffs_pol_bigQ(pol_bigQ polynomial, mpq_t *coeffs, size_t size_coeffs_array);

void set_all_coeffs_random_pol_bigQ(pol_bigQ polynomial, unsigned int max);

void set_coeff_pol_bigQ(pol_bigQ polynomial, unsigned int degree_coeff, mpq_t newValue);

void set_coeff_pol_bigQ_d(pol_bigQ polynomial, unsigned int degree_coeff, int newValue);

void set_all_coeffs_to_pol_bigQ(pol_bigQ polynomial, mpq_t value);

void set_all_coeffs_to_pol_bigQ_si(pol_bigQ polynomial, int numerator , unsigned int denominator);

void destroy_pol_bigQ(pol_bigQ polynomial);

void change_degre_pol_bigQ(pol_bigQ *polynomial, unsigned int new_degree);

void copy_pol_bigQ(pol_bigQ *res, pol_bigQ polynomial);

void print_pol_bigQ(pol_bigQ polynomial);

void add_pol_bigQ(pol_bigQ *res, pol_bigQ A, pol_bigQ B);

void sub_pol_bigQ(pol_bigQ *res, pol_bigQ A, pol_bigQ B);

void mult_pol_bigQ(pol_bigQ *res, pol_bigQ A, pol_bigQ B);

void euclideDiv_pol_bigQ(pol_bigQ *Q , pol_bigQ *R , pol_bigQ A , pol_bigQ B);



#endif //ARITHMETIC_POL_BIGQ_ARITHMETIC_H
