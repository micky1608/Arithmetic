//
// Created by root on 25/12/18.
//



#ifndef ARITHMETIC_FFT_H
#define ARITHMETIC_FFT_H

#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#include <math.h>
#include "pol_bigfloat_arithmetic.h"
#include "Complex.h"


void n_unit_roots(mpc_t X[] , unsigned int n);

void split_pol_bigfloat(pol_bigfloat *pol_even , pol_bigfloat *pol_odd , pol_bigfloat pol);
void split_pol_complex(pol_complex *pol_even , pol_complex *pol_odd , pol_complex pol);

void polynomial_naive_evaluation(mpc_t *y , pol_bigfloat P , mpc_t x);
void polynomial_naive_evaluation_C(mpc_t *y , pol_complex P , mpc_t x);


void FFT_evaluation(mpc_t *P_eval , pol_bigfloat P , mpc_t *Xi , unsigned int n);
void FFT_evaluation_C(mpc_t *P_eval , pol_complex P , mpc_t *Xi , unsigned int n);

void FFT_multiplication(mpc_t *R_eval , mpc_t *P_eval , mpc_t *Q_eval , unsigned int d);

void FFT_reverse(pol_bigfloat *R , mpc_t *R_eval , mpc_t *Xi , unsigned int n);

void mul_pol_bigfloat_FFT(pol_bigfloat *R , pol_bigfloat P , pol_bigfloat Q);

#endif //ARITHMETIC_FFT_H
