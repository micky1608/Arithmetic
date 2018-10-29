//
// Created by root on 29/10/18.
//

#ifndef ARITHMETIC_POL_BIGNUMBER_H
#define ARITHMETIC_POL_BIGNUMBER_H

#include "pol_bigfloat_arithmetic.h"
#include "pol_bigint_arithmetic.h"

void pol_bigint_to_pol_bigfloat(pol_bigfloat *res , pol_bigint A);

void pol_bigfloat_to_pol_bigint(pol_bigint *res , pol_bigfloat A);

#endif //ARITHMETIC_POL_BIGNUMBER_H
