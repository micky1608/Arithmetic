//
// Created by root on 05/10/18.
//

#ifndef MODEL_TME1_ARITHMETIC_H
#define MODEL_TME1_ARITHMETIC_H

    #include <sys/types.h>
    #include <gmp.h>


    void set32_mod_P(__uint32_t *number, __uint32_t newvalue , __uint32_t P);
    void set64_mod_P(__uint64_t *number, __uint64_t newvalue , __uint64_t P);
    void set_bigint_mod_P(mpz_t number, mpz_t newvalue, mpz_t P);

    void add32_mod_P(__uint32_t *res , __uint32_t a , __uint32_t b , __uint32_t P);
    void add64_mod_P(__uint64_t *res , __uint64_t a , __uint64_t b , __uint64_t P);
    void add_bigint_mod_P(mpz_t res , mpz_t a , mpz_t b , mpz_t P);

    void sub32_mod_P(__uint32_t *res , __uint32_t a , __uint32_t b , __uint32_t P);
    void sub64_mod_P(__uint64_t *res , __uint64_t a , __uint64_t b , __uint64_t P);
    void sub_bigint_mod_P(mpz_t res , mpz_t a , mpz_t b , mpz_t P);



#endif //MODEL_TME1_ARITHMETIC_H
