//
// Created by root on 05/10/18.
//

#include <stdio.h>
#include "arithmetic.h"

void set32_mod_P(u32 *number, u32 newvalue, u32 P) {
    *number = newvalue%P;
}

void set64_mod_P(u64 *number, u64 newvalue, u64 P) {
    *number = newvalue%P;
}

void set_bigint_mod_P(mpz_t number, mpz_t newvalue, mpz_t P) {
    mpz_mod(newvalue,newvalue,P);
    mpz_set(number , newvalue);
}

/**
 * TO USE WITH P < 2^32
 * @param res
 * @param a
 * @param b
 */
void add32_mod_P(u32 *res , u32 a , u32 b , u32 P) {
    *res = (a%P + b%P) % P;
}

/**
 * TO USE WITH P < 2^64
 * @param res
 * @param a
 * @param b
 */
void add64_mod_P(u64 *res , u64 a , u64 b , u64 P) {
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


void sub32_mod_P(u32 *res , u32 a , u32 b , u32 P) {
        if(b>a) *res = P - (b%P-a%P);
        else    *res = (a%P - b%P);
}

void sub64_mod_P(u64 *res , u64 a , u64 b , u64 P) {
    if(b>a) *res = P - (b%P-a%P);
    else    *res = (a%P - b%P);
}
void sub_bigint_mod_P(mpz_t res , mpz_t a , mpz_t b , mpz_t P) {
    mpz_sub(res,a,b);
    mpz_mod(res,res,P);
}

void mult32_mod_P(u32 *res , u32 a , u32 b , u32 P) {
    u64 temp =  ((u64)a%P)*((u64)b%P) % (u64)P;
    *res = (u32)temp;
}

void mult64_mod_P(u64 *res , u64 a , u64 b , u64 P) {
    *res = ((a%P)*(b%P))%P;
}

void mult_bigint_mod_P(mpz_t res , mpz_t a , mpz_t b , mpz_t P) {
    mpz_mul(res,a,b);
    mpz_mod(res,res,P);
}

void div32_mod_P(u32 *res , u32 a , u32 b , u32 P) {
    *res = (a%P)/(b%P);
}

void div64_mod_P(u64 *res , u64 a , u64 b , u64 P) {
    *res = (a%P)/(b%P);
}

void div_bigint_mod_P(mpz_t res , mpz_t a , mpz_t b , mpz_t P) {
    mpz_div(res,a,b);
    mpz_mod(res,res,P);
}





