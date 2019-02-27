//
// Created by root on 13/10/18.
//

#ifndef ARITHMETIC_UTIL_H
#define ARITHMETIC_UTIL_H

#include <math.h>
#include <gmp.h>

typedef enum bool { TRUE, FALSE } bool;

int max(int a , int b);

double log_base_2(double x);

void mpq_mul_si(mpq_t res , mpq_t x , int lambda_num , unsigned int lambda_den);

void mpq_pow(mpq_t res , mpq_t x , int p);

void swap(double *a , double *b);

void swap_long(long *a , long *b);

#endif //ARITHMETIC_UTIL_H
