//
// Created by root on 05/10/18.
//

#ifndef MODEL_TME1_ARITHMETIC_H
#define MODEL_TME1_ARITHMETIC_H

#include <sys/types.h>
#include <gmp.h>

    typedef __uint32_t u32;
    typedef __uint64_t u64;


    void set32_mod_P(u32 *number, u32 newvalue , u32 P);
    void set64_mod_P(u64 *number, u64 newvalue , u64 P);
    void set_bigint_mod_P(mpz_t number, mpz_t newvalue, mpz_t P);

    void add32_mod_P(u32 *res , u32 a , u32 b , u32 P);
    void add64_mod_P(u64 *res , u64 a , u64 b , u64 P);
    void add_bigint_mod_P(mpz_t res , mpz_t a , mpz_t b , mpz_t P);

    void sub32_mod_P(u32 *res , u32 a , u32 b , u32 P);
    void sub64_mod_P(u64 *res , u64 a , u64 b , u64 P);
    void sub_bigint_mod_P(mpz_t res , mpz_t a , mpz_t b , mpz_t P);

    void div64_mod_P(int *res , int a , int b , int P);
    void div_bigint_mod_P(mpz_t res , mpz_t a , mpz_t b , mpz_t P);

    void EEA64(int *U , int *V , int A , int B);
    void EEA_bigint(mpz_t U , mpz_t V , mpz_t A , mpz_t B);

    void inv64_mod_P(int *inv , int A , int P);
    void inv_bigint_mod_P(mpz_t inv , mpz_t A , mpz_t P);



#endif //MODEL_TME1_ARITHMETIC_H
