//
// Created by root on 24/12/18.
//
/*
A complex rounding mode is of the form MPC_RNDxy where x and y are one of :
 - N (to nearest)
 - Z (towards zero)
 - U (towards plus infinity)
 - D (towards minus infinity).
 The first letter refers to the rounding mode for the real part, and the second one for the imaginary part.
 For example MPC_RNDZU indicates to round the real part towards zero, and the imaginary part towards plus infinity.
*/


#ifndef ARITHMETIC_COMPLEX_H
#define ARITHMETIC_COMPLEX_H

#include <mpc.h>
#include <mpfr.h>
#include <stdio.h>
#include <malloc.h>

typedef struct pol_complex {
    unsigned int degree;
    mpc_t *coeffs;
} pol_complex;

void print_complex(mpc_t z , char *name);

void init_pol_complex(pol_complex *pol , unsigned int degree);

void set_coeffs_pol_complex(pol_complex *pol , mpc_t *new_coeffs , unsigned int nb_coeffs);

void change_degre_pol_complex(pol_complex *pol, unsigned int new_degree);

void destroy_pol_complex(pol_complex pol);



#endif //ARITHMETIC_COMPLEX_H
