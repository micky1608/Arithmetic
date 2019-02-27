//
// Created by root on 13/10/18.
//

#include "Util.h"


int max(int a , int b) {
    if(a >= b) return a;
    return b;
}

/* ********************************************************************************************************************** */

double log_base_2(double x) {
    return (log(x) / log(2));
}

/* ********************************************************************************************************************** */

void swap(double *a , double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

/* ********************************************************************************************************************** */

void swap_long(long *a , long *b) {
    long temp = *a;
    *a = *b;
    *b = temp;
}

/* ********************************************************************************************************************** */

void mpq_mul_si(mpq_t res , mpq_t x , int lambda_num , unsigned int lambda_den) {
    mpz_t numerator, denominator,gcd;

    mpz_inits(numerator , denominator , gcd , NULL);

    mpq_get_num(numerator , x);
    mpq_get_den(denominator , x);

    mpz_mul_si(numerator , numerator , lambda_num);
    mpz_mul_si(denominator , denominator , lambda_den);

    mpz_gcd(gcd , numerator , denominator);
    mpz_div(numerator , numerator , gcd);
    mpz_div(denominator , denominator , gcd);

    mpq_set_num(res , numerator);
    mpq_set_den(res , denominator);

    mpz_clears(numerator , denominator , gcd , NULL);
}

/* ********************************************************************************************************************** */

void mpq_pow(mpq_t res , mpq_t x , int p) {

    mpq_set(res , x);
    for(int i=1 ; i<p ; i++)
        mpq_mul(res , res , x);


}
