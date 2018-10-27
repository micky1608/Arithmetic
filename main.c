//
// Created by root on 05/10/18.
//


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include "arithmetic.h"
#include "polynomial_arithmetic.h"

int main () {

    srand(time(NULL));

    mpz_t e,f,PM,resM,P;
    mpz_inits(resM,P,(mpz_t *)NULL);

    uint32_t P32 = 2147483647;
    uint64_t P64 = 2305843009213693951ULL;
    mpz_init_set_str(PM, "618970019642690137449562111",10);

    int res;

    uint32_t a = P32 - 1;
    uint32_t b = 10;
    uint32_t res32;

    uint64_t c = P64 -1 ;
    uint64_t d = 10;
    uint64_t res64;

    mpz_init_set_str(e , "2000000000000000000000000", 10);
    mpz_init_set_str(f , "2000000000000000000000000", 10);


    add32_mod_P(&res32 , a , b , P32);
    printf("Result in Z/P32.Z :\ta + b = %u\n" , res32);

    add64_mod_P(&res64 , c ,d , P64);
    printf("Result in Z/P64.Z :\tc + d = %lu\n" , res64);

    add_bigint_mod_P(resM,e,f,PM);

    printf("Result in Z/PM.Z :\te + f = ");
    mpz_out_str(stdout,10,resM);

    printf("\n\n*******************************************\n\n");

    set32_mod_P(&a , (P32+15) , P32);
    printf("Value in Z/P32.Z of P32+15 : %u\n",a);

    set64_mod_P(&c , (P64+3) , P64);
    printf("Value in Z/P64.Z of P64+3 : %lu\n",c);

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
    printf("Value in Z/P64.Z : c = %lu\n",c);
    printf("Value in Z/P64.Z : d = %lu\n",d);
    printf("Result in Z/P64.Z :\tc - d = %lu " , res64);
    if(res64 == P64-1)
        printf("= P64 - 1\n");
    else
        printf("!= P64 - 1\n");




    printf("\n\n*******************************************\n\n");

    div64_mod_P(&res , 2 , 4 , 5);
    printf("2/4 mod 5 = %d\n" , res);

    div64_mod_P(&res , 4 , 2 , 5);
    printf("4/2 mod 5 = %d\n" , res);

    div64_mod_P(&res , 1 , 3 , 5);
    printf("2/4 mod 5 = %d\n" , res);


    mpz_set_d(e,2);
    mpz_set_d(f,4);
    mpz_set_d(P,5);

    inv_bigint_mod_P(resM , e , P);
    printf("Inv of 2 mod 5 : ");
    mpz_out_str(stdout,10,resM);
    printf("\n");

    inv_bigint_mod_P(resM , f , P);
    printf("Inv of 4 mod 5 : ");
    mpz_out_str(stdout,10,resM);
    printf("\n");


    printf("\n\n*******************************************\n\n");


    Pol_M my_polynomial;
    init_Pol_M(&my_polynomial , 5);
    set_all_coeffs_random_Pol_M(my_polynomial , 100);

    printf("Random Polynomial of degree 5 with coefficients between -100 and 100 : \t");
    print_Pol_m(my_polynomial);

    mpz_set_d(e,-20);
    set_coeff_Pol_M(my_polynomial,3,e);
    printf("Same polynomial with the coefficient od degree 3 = -20 : \t\t\t\t");
    print_Pol_m(my_polynomial);

    change_degre_Pol_M(&my_polynomial , 6);
    printf("Same polynomial extended to the degree 6 : \t\t\t\t\t\t\t\t");
    print_Pol_m(my_polynomial);

    change_degre_Pol_M(&my_polynomial , 4);
    printf("Same polynomial shrank to the degree 4 : \t\t\t\t\t\t\t\t");
    print_Pol_m(my_polynomial);


    printf("\n\n*******************************************\n\n");

    Pol_M other_polynomial , res_polynomial;
    init_Pol_M(&other_polynomial , 5);
    init_Pol_M(&res_polynomial , 2);

    set_all_coeffs_random_Pol_M(my_polynomial , 100);
    set_all_coeffs_random_Pol_M(other_polynomial , 100);

    add_Pol_M(&res_polynomial , other_polynomial , my_polynomial);

    printf("A : \t");
    print_Pol_m(my_polynomial);

    printf("B : \t");
    print_Pol_m(other_polynomial);

    printf("A + B : ");
    print_Pol_m(res_polynomial);

    printf("\n");

    change_degre_Pol_M(&my_polynomial , 2);
    change_degre_Pol_M(&other_polynomial , 3);

    set_all_coeffs_random_Pol_M(my_polynomial , 5);
    set_all_coeffs_random_Pol_M(other_polynomial , 5);

    mult_Pol_M(&res_polynomial , other_polynomial , my_polynomial);

    printf("A : \t");
    print_Pol_m(my_polynomial);

    printf("B : \t");
    print_Pol_m(other_polynomial);

    printf("A*B : \t");
    print_Pol_m(res_polynomial);


    printf("\n\n*******************************************\n\n");

    Pol_M remainder_pol;
    init_Pol_M(&remainder_pol , 0);

//    change_degre_Pol_M(&my_polynomial , 5);
//    change_degre_Pol_M(&other_polynomial , 2);
//
//    set_all_coeffs_random_Pol_M(my_polynomial , 10);
//    set_all_coeffs_random_Pol_M(other_polynomial , 10);
//
//    printf("\tEuclidean division\n");
//    printf("A : \t");
//    print_Pol_m(my_polynomial);
//
//    printf("B : \t");
//    print_Pol_m(other_polynomial);
//
//    euclide_div_Pol_M(&res_polynomial , &remainder_pol , my_polynomial , other_polynomial);
//
//
//
//    printf("\n\tResult\n");
//    printf("Q : \t");
//    print_Pol_m(res_polynomial);
//
//    printf("R : \t");
//    print_Pol_m(remainder_pol);


    destroy_Pol_M(my_polynomial);
    destroy_Pol_M(other_polynomial);
    destroy_Pol_M(res_polynomial);
    destroy_Pol_M(remainder_pol);

    mpz_clears(e,f,PM,resM,P,(mpz_t *)NULL);

    return 0;
}
