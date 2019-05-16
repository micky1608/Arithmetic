//
// Created by micky on 27/02/19.
//

#ifndef ARITHMETIC_TEST_H
#define ARITHMETIC_TEST_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <mpc.h>

#include "arithmetic.h"
#include "pol_bigint_arithmetic.h"
#include "pol_bigfloat_arithmetic.h"
#include "pol_bigQ_arithmetic.h"
#include "Euclidean.h"
#include "matrix_bigint_arithmetic.h"
#include "matrix_bigQ_arithmetic.h"
#include "matrix_double.h"
#include "CRT_bigint.h"
#include "Lagrange.h"
#include "Complex.h"
#include "FFT.h"
#include "pol.h"
#include "matrix_pol.h"
#include "matrix_pol_double.h"
#include "pol_double.h"
#include "tree.h"
#include "pol_factorization.h"
#include "serie.h"

void test_pol();
void test_halfGCD();
void test_horner();
void test_matrix_pol();
void test_fastEuclide();
void test_pol_fact();
void test_serie();


#endif //ARITHMETIC_TEST_H
