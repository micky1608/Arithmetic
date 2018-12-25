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

void print_complex(mpc_t z , char *name);

#endif //ARITHMETIC_COMPLEX_H
