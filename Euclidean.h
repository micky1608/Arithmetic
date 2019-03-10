//
// Created by root on 29/10/18.
//

#ifndef ARITHMETIC_EUCLIDEAN_H
#define ARITHMETIC_EUCLIDEAN_H

#include <math.h>

#include "pol_bigint_arithmetic.h"
#include "pol_bigfloat_arithmetic.h"
#include "pol_bignumber.h"
#include "pol.h"
#include "matrix.h"
#include "matrix_pol.h"

void euclideDiv_pol_bignumber(pol_bigfloat *Q , pol_bigfloat *R , pol_bigint A , pol_bigint B);

void halfGCD(matrix_pol *Mgcd , pol A , pol B);



#endif //ARITHMETIC_EUCLIDEAN_H
