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


void split_pol_bigfloat(pol_bigfloat *pol_even , pol_bigfloat *pol_odd , pol_bigfloat pol);

void n_unit_roots(mpc_t X[] , unsigned int n);

void FFT_evaluation(mpc_t *P_eval , pol_bigfloat P , mpc_t omega , unsigned int n);

void FFT_multiplication(mpc_t *R_eval , mpc_t *P_eval , mpc_t Q_eval);

void FFT_reverse(pol_bigfloat *R , mpc_t *R_eval , mpc_t omega , unsigned int n);

void mul_pol_bigfloat_FFT(pol_bigfloat *R , pol_bigfloat P , pol_bigfloat Q);

#endif //ARITHMETIC_FFT_H
