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
    mpz_t U,V,GCD,ONE, A_MOD_P, INV_B_MOD_P;
    mpz_inits(U,V,GCD,ONE,A_MOD_P,INV_B_MOD_P,0);
    mpz_set_d(ONE,1);
    mpz_mod(A_MOD_P , a , P);

    EEA_bigint_mod_P(&U , &V , &GCD , ONE , b);

    if(!mpz_cmp_d(GCD,1)) {

        mpz_set(INV_B_MOD_P , U);
        mpz_mod(INV_B_MOD_P , INV_B_MOD_P , P);
        mpz_add(INV_B_MOD_P , INV_B_MOD_P , P);
        mpz_mod(INV_B_MOD_P , INV_B_MOD_P , P);

        mpz_mul(res , A_MOD_P, INV_B_MOD_P);
    }

    mpz_clears(U,V,GCD,ONE,A_MOD_P,INV_B_MOD_P,0);

}


void EEA_bigint_mod_P(mpz_t *U , mpz_t *V , mpz_t *GCD , mpz_t A , mpz_t B) {

    mpz_t RI , RI_1; // GCD is R(i-1)
    mpz_t UI , UI_1; // U is U(i-1)
    mpz_t VI , VI_1; // V is V(i-1)
    mpz_t Q;
    int i=1;
    mpz_inits(RI,RI_1,UI,UI_1,VI,VI_1,Q,0);

    mpz_set(*GCD , A);
    mpz_set(RI , B);
    mpz_set_d(*U , 1);
    mpz_set_d(UI , 0);
    mpz_set_d(*V , 0);
    mpz_set_d(VI , 1);

    while(mpz_cmp_d(RI , 0)) {
        mpz_div(Q , *GCD , RI);
        mpz_mod(RI_1 , *GCD , RI);

        mpz_set(UI_1 , Q);
        mpz_neg(UI_1 , UI_1);
        mpz_mul(UI_1 , UI_1 , UI);
        mpz_add(UI_1 , UI_1 , *U);

        mpz_set(VI_1 , Q);
        mpz_neg(VI_1 , VI_1);
        mpz_mul(VI_1 , VI_1 , VI);
        mpz_add(VI_1 , VI_1 , *V);

        i++;
        mpz_set(*GCD , RI);
        mpz_set(RI , RI_1);

        mpz_set(*GCD , RI);
        mpz_set(RI , RI_1);

        mpz_set(*U , UI);
        mpz_set(UI , UI_1);

        mpz_set(*V , VI);
        mpz_set(VI , VI_1);

    }

    mpz_clears(RI,RI_1,UI,UI_1,VI,VI_1,Q,0);

}






