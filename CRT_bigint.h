//
// Created by root on 23/11/18.
//



#ifndef ARITHMETIC_CRT_BIGINT_H
#define ARITHMETIC_CRT_BIGINT_H

#include <gmp.h>
#include "arithmetic.h"
#include "Util.h"


typedef struct crt {
    mpz_t Ai;
    mpz_t Ni;
} crt_bigint;


void init_CRT_bigint(crt_bigint *crt , mpz_t A , mpz_t N);

void init_CRT_bigint_d(crt_bigint *crt , int A , int N);

void destroy(crt_bigint crt);

void print_CRT_bigint(crt_bigint CRT);

void print_CRT_bigint_arg(crt_bigint *arg , unsigned int nb_arg);

void solve_CRT_bigint(crt_bigint *res , crt_bigint *arg , unsigned int nb_arg);

void solve_CRT_bigint_DAC(crt_bigint *res , crt_bigint *arg , unsigned int nb_arg);



#endif //ARITHMETIC_CRT_BIGINT_H
