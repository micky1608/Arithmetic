//
// Created by root on 12/10/18.
//

#ifndef ARITHMETIC_POLYNOMIAL_ARITHMETIC_H
#define ARITHMETIC_POLYNOMIAL_ARITHMETIC_H

#include <gmp.h>
#include "Util.h"

typedef struct Pol_M {
    unsigned int degree;
    mpz_t *coeffs;
} Pol_M;


void init_Pol_M(Pol_M *polynomial , unsigned int degree);

void set_all_coeffs_Pol_M(Pol_M polynomial , mpz_t *coeffs , size_t size_coeffs_array);

void set_all_coeffs_random_Pol_M(Pol_M polynomial , unsigned int max);

void set_coeff_Pol_M(Pol_M polynomial , unsigned int degree_coeff , mpz_t newValue);

void set_all_coeffs_to_Pol_M(Pol_M polynomial , mpz_t value);

void destroy_Pol_M(Pol_M polynomial);

void change_degre_Pol_M(Pol_M *polynomial , unsigned int new_degree);

void copy_Pol_M(Pol_M *res , Pol_M polynomial);

void print_Pol_m(Pol_M polynomial);

void add_Pol_M(Pol_M *res , Pol_M A , Pol_M B);

void mult_Pol_M(Pol_M *res , Pol_M A , Pol_M B);

void euclide_div_Pol_M_Mod_P(Pol_M *Q , Pol_M *R , Pol_M A , Pol_M B , mpz_t P);



#endif //ARITHMETIC_POLYNOMIAL_ARITHMETIC_H
