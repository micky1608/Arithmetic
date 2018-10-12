//
// Created by root on 12/10/18.
//

#ifndef ARITHMETIC_POLYNOMIAL_ARITHMETIC_H
#define ARITHMETIC_POLYNOMIAL_ARITHMETIC_H

#include <gmp.h>

typedef struct Pol_M {
    unsigned int degree;
    mpz_t *coeffs;
} Pol_M;


void init_Pol_M(Pol_M *polynomial , unsigned int degree);

void set_all_coeffs_Pol_M(Pol_M polynomial , mpz_t *coeffs);

void set_all_coeffs_random_Pol_M(Pol_M polynomial , unsigned int max);

void set_coeff_Pol_M(Pol_M polynomial , unsigned int degree_coeff , mpz_t newValue);

void destroy_Pol_M(Pol_M polynomial);

void change_degre_Pol_M(Pol_M *polynomial , unsigned int new_degree);

void print_Pol_m(Pol_M polynomial);



#endif //ARITHMETIC_POLYNOMIAL_ARITHMETIC_H
