//
// Created by root on 30/11/18.
//

#ifndef ARITHMETIC_LAGRANGE_H
#define ARITHMETIC_LAGRANGE_H

#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include "pol_bigQ_arithmetic.h"

typedef struct pol_bigQ_value {
    mpq_t Xi;
    mpq_t Yi;
} pol_bigQ_value;

void init_polbigQ_value(pol_bigQ_value *polBigQValue , mpq_t Xi , mpq_t Yi);

void init_polbigQ_value_si(pol_bigQ_value *polBigQValue , int Xi_num , unsigned int Xi_den, int Yi_num , unsigned int Yi_den);

void destroy_polbigQ_value(pol_bigQ_value polBigQValue);

void lagrange_interpolation(pol_bigQ *res , pol_bigQ_value *values , unsigned int nb_values);

void print_polbigQ_value(pol_bigQ_value polValue);

void print_lagrange_constaints(pol_bigQ_value *constraints , unsigned int nb_constraint);

#endif //ARITHMETIC_LAGRANGE_H
