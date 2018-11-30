//
// Created by root on 30/11/18.
//

#ifndef ARITHMETIC_LAGRANGE_H
#define ARITHMETIC_LAGRANGE_H

#include <gmp.h>
#include "pol_bigQ_arithmetic.h"

typedef struct pol_bigQ_value {
    mpq_t Xi;
    mpq_t Yi;
} pol_value;

void lagrange_interpolation(pol_bigQ *res , pol_value *values , unsigned int nb_values);

#endif //ARITHMETIC_LAGRANGE_H
