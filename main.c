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


    uint32_t a = P32 - 1;
    uint32_t b = 10;
    uint32_t res32;

    uint64_t c = P64 -1 ;
    uint64_t d = 10;
    uint64_t res64;

    mpz_init_set_str(e , "2000000000000000000000000", 10);
    mpz_init_set_str(f , "2000000000000000000000000", 10);


    add32_mod_P(&res32 , a , b , P32);
    printf("Result in Z/P32.Z :\ta + b = %d\n" , res32);

    add64_mod_P(&res64 , c ,d , P64);
    printf("Result in Z/P64.Z :\tc + d = %lu\n" , res64);

    add_bigint_mod_P(resM,e,f,PM);

    printf("Result in Z/PM.Z :\te + f = ");
    mpz_out_str(stdout,10,resM);

    printf("\n\n*******************************************\n\n");

    set32_mod_P(&a , (P32+15) , P32);
    printf("Value in Z/P32.Z of P32+15 : %d\n",a);

    set64_mod_P(&c , (P64+3) , P64);
    printf("Value in Z/P64.Z of P64+3 : %d\n",c);

    set_bigint_mod_P(e,PM,PM);
    printf("Value in Z/PM.Z of PM : ");
    mpz_out_str(stdout,10,e);
    printf("\n");



    mpz_clears(e,f,PM,resM,0);


    return 0;
}
