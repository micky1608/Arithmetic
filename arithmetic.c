//
// Created by root on 05/10/18.
//

#include "arithmetic.h"

void set32_mod_P(__uint32_t *number, __uint32_t newvalue, __uint32_t P) {

}

void set64_mod_P(__uint64_t *number, __uint64_t newvalue, __uint64_t P) {

}

void set_bigint_mod_P(mpz_t number, mpz_t newvalue, mpz_t P) {

}

/**
 * TO USE WITH P < 2^32
 * @param res
 * @param a
 * @param b
 */
void add32_mod_P(__uint32_t *res , __uint32_t a , __uint32_t b , __uint32_t P) {
    *res = (a%P + b%P) % P;
}

/**
 * TO USE WITH P < 2^64
 * @param res
 * @param a
 * @param b
 */
void add64_mod_P(__uint64_t *res , __uint64_t a , __uint64_t b , __uint64_t P) {
    *res = (a%P + b%P) % P;
}

/**
 * TO USE WITH P >= 2^64
 * @param res
 * @param a
 * @param b
 */
void add_bigint_mod_P(mpz_t res , mpz_t a , mpz_t b , mpz_t P) {
    mpz_add(res,a,b);
    mpz_mod(res,res,P);
}



