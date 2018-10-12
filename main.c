//
// Created by root on 05/10/18.
//


#include <stdio.h>
#include <stdint.h>
#include "arithmetic.h"
#include "polynomial_arithmetic.h"

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

    printf("\n\n*******************************************\n\n");

    set32_mod_P(&a, 0 , P32);
    set32_mod_P(&b, 1 , P32);
    sub32_mod_P(&res32,a,b,P32);
    printf("Value in Z/P32.Z : a = %u\n",a);
    printf("Value in Z/P32.Z : b = %u\n",b);
    printf("Result in Z/P32.Z :\ta - b = %u " , res32);
    if(res32 == P32-1)
        printf("= P32 - 1\n");
    else
        printf("!= P32 - 1\n");

    printf("\n\n*******************************************\n\n");

    set64_mod_P(&c, 0 , P64);
    set64_mod_P(&d, 1 , P64);
    sub64_mod_P(&res64,c,d,P64);
    printf("Value in Z/P64.Z : c = %u\n",c);
    printf("Value in Z/P64.Z : d = %u\n",d);
    printf("Result in Z/P64.Z :\tc - d = %u " , res64);
    if(res64 == P64-1)
        printf("= P64 - 1\n");
    else
        printf("!= P64 - 1\n");


    mpz_clears(e,f,PM,resM,0);



    /************************************************************************************************************************************** */

    printf("\n\n*******************************************\n\n");

    Pol_M my_polynomial;
    init_Pol_M(&my_polynomial , 5);
    set_all_coeffs_random_Pol_M(my_polynomial , 100);
    print_Pol_m(my_polynomial);



    return 0;
}
