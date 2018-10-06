//
// Created by root on 05/10/18.
//


#include <stdio.h>
#include <stdint.h>
#include "arithmetic.h"

int main () {

    mpz_t e,f,PM,resM;
    mpz_init(resM);

    uint32_t P32 = 2147483647;
    uint64_t P64 = 2305843009213693951ULL;
    mpz_init_set_str(PM, "618970019642690137449562111",10);


    uint32_t a = 10000;
    uint32_t b = 5000;
    uint32_t res32;

    uint64_t c = P64 * 5;
    uint64_t d = P64 * 10;
    uint64_t res64;

    mpz_init_set_str(e , "2000000000000000000000000", 10);
    mpz_init_set_str(f , "2000000000000000000000000", 10);


    add32_mod_P(&res32 , a , b , P32);
    printf("Result 32: %d\n" , res32);

    add64_mod_P(&res64 , c ,d , P64);
    printf("Result 64: %lu\n" , res64);

    add_bigint_mod_P(resM,e,f,PM);

    printf("Result M : ");
    mpz_out_str(stdout,10,resM);



    mpz_clears(e,f,PM,0);


    return 0;
}
