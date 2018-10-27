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



void div64_mod_P(int *res , int A , int B , int P) {
    if(B == 0) {
        perror("Divide by 0 !");
        return;
    }

    int inv_B_mod_P;
    inv64_mod_P(&inv_B_mod_P , B , P);
    *res = ((A%P) * inv_B_mod_P) % P;
}



void div_bigint_mod_P(mpz_t res , mpz_t a , mpz_t b , mpz_t P) {
    mpz_t inv_B_mod_P;
    mpz_init(inv_B_mod_P);

    inv_bigint_mod_P(inv_B_mod_P , b , P);
    mpz_mod(res , a , P);
    mpz_mul(res , res , inv_B_mod_P);
    mpz_mod(res , res , P);

    mpz_clear(inv_B_mod_P);
}

void EEA64(int *U , int *V , int A , int B) {
    int RI_previous , RI , RI_next;
    int UI , UI_next; // U is U(i-1)
    int VI , VI_next; // V is V(i-1)
    int Q;

    RI_previous = A;
    RI = B;
    *U = 1;
    UI = 0;
    *V = 0;
    VI = 1;

    while (RI != 0) {
        Q = RI_previous / RI;
        RI_next = RI_previous - Q*RI;
        UI_next = *U - Q*UI;
        VI_next = *V - Q*VI;

        RI_previous = RI;
        RI = RI_next;
        *U = UI;
        UI = UI_next;
        *V = VI;
        VI = VI_next;
    }
}


void EEA_bigint(mpz_t U , mpz_t V , mpz_t A , mpz_t B) {

    mpz_t RI_previous , RI , RI_next;
    mpz_t UI, UI_next; // U is U(i-1)
    mpz_t VI, VI_next; // V is V(i-1)
    mpz_t Q;

    //mpz_inits(RI_previous , RI , RI_next , UI , UI_next , VI , VI_next , Q , NULL);

    mpz_init(RI_previous);
    mpz_init(RI);
    mpz_init(RI_next);
    mpz_init(UI);
    mpz_init(UI_next);
    mpz_init(VI);
    mpz_init(VI_next);
    mpz_init(Q);


    mpz_set(RI_previous , A);
    mpz_set(RI , B);
    mpz_set_d(U , 1);
    mpz_set_d(UI , 0);
    mpz_set_d(V , 0);
    mpz_set_d(VI , 1);

    while(mpz_cmp_d(RI , 0)) {

        mpz_div(Q , RI_previous , RI);

        mpz_set(RI_next , Q);
        mpz_neg(RI_next , RI_next);
        mpz_mul(RI_next , RI_next , RI);
        mpz_add(RI_next , RI_next , RI_previous);

        mpz_set(UI_next , Q);
        mpz_neg(UI_next , UI_next);
        mpz_mul(UI_next , UI_next , UI);
        mpz_add(UI_next , UI_next , U);

        mpz_set(VI_next , Q);
        mpz_neg(VI_next , VI_next);
        mpz_mul(VI_next , VI_next , VI);
        mpz_add(VI_next , VI_next , V);


        mpz_set(RI_previous , RI);
        mpz_set(RI , RI_next);

        mpz_set(U , UI);
        mpz_set(UI , UI_next);

        mpz_set(V , VI);
        mpz_set(VI , VI_next);
    }

    mpz_clears(RI_previous , RI , RI_next , UI , UI_next , VI , VI_next , Q , NULL);

}

void inv64_mod_P(int *inv , int A , int P) {
    int U,V;
    EEA64(&U , &V , A , P);
    *inv = (U%P + P) % P;
}

void inv_bigint_mod_P(mpz_t inv , mpz_t A , mpz_t P) {
    mpz_t U,V;
    mpz_inits(U,V,NULL);

    EEA_bigint(U,V,A,P);

    mpz_mod(inv , U , P);
    mpz_add(inv , inv , P);
    mpz_mod(inv , inv , P);

    mpz_clears(U,V,NULL);
}










