//
// Created by root on 05/10/18.
//


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include "arithmetic.h"
#include "pol_bigint_arithmetic.h"
#include "pol_bigfloat_arithmetic.h"
#include "pol_bigQ_arithmetic.h"
#include "Euclidean.h"
#include "matrix_bigint_arithmetic.h"
#include "matrix_bigQ_arithmetic.h"
#include "matrix_double.h"
#include "CRT_bigint.h"
#include "Lagrange.h"

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


    pol_bigint my_polynomial;
    init_pol_bigint(&my_polynomial, 5);
    set_all_coeffs_random_pol_bigint(my_polynomial, 100);

    printf("Random Polynomial of degree 5 with coefficients between -100 and 100 : \t");
    print_pol_bigint(my_polynomial);

    mpz_set_d(e,-20);
    set_coeff_pol_bigint(my_polynomial, 3, e);
    printf("Same polynomial with the coefficient od degree 3 = -20 : \t\t\t\t");
    print_pol_bigint(my_polynomial);

    change_degre_pol_bigint(&my_polynomial, 6);
    printf("Same polynomial extended to the degree 6 : \t\t\t\t\t\t\t\t");
    print_pol_bigint(my_polynomial);

    change_degre_pol_bigint(&my_polynomial, 4);
    printf("Same polynomial shrank to the degree 4 : \t\t\t\t\t\t\t\t");
    print_pol_bigint(my_polynomial);


    printf("\n\n*******************************************\n\n");

    pol_bigint other_polynomial , res_polynomial;
    init_pol_bigint(&other_polynomial, 5);
    init_pol_bigint(&res_polynomial, 2);

    set_all_coeffs_random_pol_bigint(my_polynomial, 100);
    set_all_coeffs_random_pol_bigint(other_polynomial, 100);

    add_pol_bigint(&res_polynomial, other_polynomial, my_polynomial);

    printf("A : \t");
    print_pol_bigint(my_polynomial);

    printf("B : \t");
    print_pol_bigint(other_polynomial);

    printf("A + B : ");
    print_pol_bigint(res_polynomial);

    sub_pol_bigint(&res_polynomial, my_polynomial, other_polynomial);
    printf("A - B : ");
    print_pol_bigint(res_polynomial);

    printf("\n");

    change_degre_pol_bigint(&my_polynomial, 2);
    change_degre_pol_bigint(&other_polynomial, 3);

    set_all_coeffs_random_pol_bigint(my_polynomial, 5);
    set_all_coeffs_random_pol_bigint(other_polynomial, 5);

    mult_pol_bigint(&res_polynomial, other_polynomial, my_polynomial);

    printf("A : \t");
    print_pol_bigint(my_polynomial);

    printf("B : \t");
    print_pol_bigint(other_polynomial);

    printf("A*B : \t");
    print_pol_bigint(res_polynomial);


    printf("\n\n*******************************************\n\n");

    pol_bigfloat quotient_pol_bigfloat, remainder_pol_bigfloat;

    init_pol_bigfloat(&quotient_pol_bigfloat, 0);
    init_pol_bigfloat(&remainder_pol_bigfloat, 0);

    change_degre_pol_bigint(&my_polynomial , 2);
    change_degre_pol_bigint(&other_polynomial , 1);

    set_coeff_pol_bigint_d(my_polynomial , 0 , 3);
    set_coeff_pol_bigint_d(my_polynomial , 1 , 2);
    set_coeff_pol_bigint_d(my_polynomial , 2 , 1);

    set_coeff_pol_bigint_d(other_polynomial , 0 , -4);
    set_coeff_pol_bigint_d(other_polynomial , 1 , 1);


    printf("\t\t***** Euclidean division *****\n\n");
    printf("A : \t");
    print_pol_bigint(my_polynomial);

    printf("B : \t");
    print_pol_bigint(other_polynomial);

    euclideDiv_pol_bignumber(&quotient_pol_bigfloat , &remainder_pol_bigfloat , my_polynomial , other_polynomial);

    printf("Q : \t");
    print_pol_bigfloat(quotient_pol_bigfloat);

    printf("R : \t");
    print_pol_bigfloat(remainder_pol_bigfloat);

    change_degre_pol_bigint(&my_polynomial , 5);
    change_degre_pol_bigint(&other_polynomial , 2);

    set_coeff_pol_bigint_d(my_polynomial , 0 , 0);
    set_coeff_pol_bigint_d(my_polynomial , 1 , -2);
    set_coeff_pol_bigint_d(my_polynomial , 2 , 3);
    set_coeff_pol_bigint_d(my_polynomial , 3 , -1);
    set_coeff_pol_bigint_d(my_polynomial , 4 , -1);
    set_coeff_pol_bigint_d(my_polynomial , 5 , 1);

    set_coeff_pol_bigint_d(other_polynomial , 0 , 1);
    set_coeff_pol_bigint_d(other_polynomial , 1 , -1);
    set_coeff_pol_bigint_d(other_polynomial , 2 , 1);

    destroy_pol_bigfloat(quotient_pol_bigfloat);
    destroy_pol_bigfloat(remainder_pol_bigfloat);
    init_pol_bigfloat(&quotient_pol_bigfloat , 0);
    init_pol_bigfloat(&remainder_pol_bigfloat , 0);

    printf("\n");
    printf("A : \t");
    print_pol_bigint(my_polynomial);

    printf("B : \t");
    print_pol_bigint(other_polynomial);

    euclideDiv_pol_bignumber(&quotient_pol_bigfloat , &remainder_pol_bigfloat , my_polynomial , other_polynomial);

    printf("Q : \t");
    print_pol_bigfloat(quotient_pol_bigfloat);

    printf("R : \t");
    print_pol_bigfloat(remainder_pol_bigfloat);


    printf("\n\n*******************************************\n\n");

    printf("\t\t***** Polynomial multiplication *****\n\n");

    pol_bigint F , G;
    init_pol_bigint(&F, 3);
    init_pol_bigint(&G, 3);

    mpz_set_d(e , 1);
    mpz_set_d(f , 2);
    set_all_coeffs_to_pol_bigint(F, e);
    set_all_coeffs_to_pol_bigint(G, f);

    destroy_pol_bigint(res_polynomial);
    init_pol_bigint(&res_polynomial, 6);
    karatsuba_pol_bigint(&res_polynomial, F, G);

    printf("F : ");
    print_pol_bigint(F);
    printf("G : ");
    print_pol_bigint(G);

    printf("F * G (karatsuba): \t");
    print_pol_bigint(res_polynomial);

    mult_pol_bigint(&res_polynomial, F, G);
    printf("F * G (naive)\t: \t");
    print_pol_bigint(res_polynomial);

    printf("\n");

    change_degre_pol_bigint(&F, 7);
    change_degre_pol_bigint(&G, 7);
    set_all_coeffs_to_pol_bigint(F, e);
    set_all_coeffs_to_pol_bigint(G, f);

    destroy_pol_bigint(res_polynomial);
    init_pol_bigint(&res_polynomial, 14);
    karatsuba_pol_bigint(&res_polynomial, F, G);

    printf("F : ");
    print_pol_bigint(F);
    printf("G : ");
    print_pol_bigint(G);

    printf("F * G (karatsuba): \t");
    print_pol_bigint(res_polynomial);

    mult_pol_bigint(&res_polynomial, F, G);
    printf("F * G (naive)\t: \t");
    print_pol_bigint(res_polynomial);


    destroy_pol_bigint(my_polynomial);
    destroy_pol_bigint(other_polynomial);
    destroy_pol_bigint(res_polynomial);
    destroy_pol_bigfloat(remainder_pol_bigfloat);

    printf("\n\n*******************************************\n\n");

    matrix_bigint matrix;

    init_matrix_bigint(&matrix , 3 , 4);
    set_coeff_matrix_bigint_d(matrix , 1 , 3 , 2);
    print_matrix_bigint(matrix);

    set_all_coeffs_random_matrix_bigint(matrix , 100);
    print_matrix_bigint(matrix);

    change_nb_line_matrix_bigint(&matrix , 6);
    print_matrix_bigint(matrix);

    change_nb_line_matrix_bigint(&matrix , 2);
    print_matrix_bigint(matrix);

    change_nb_col_matrix_bigint(&matrix , 6);
    print_matrix_bigint(matrix);

    change_nb_col_matrix_bigint(&matrix , 2);
    print_matrix_bigint(matrix);

    change_dim_matrix_bigint(&matrix , 3 , 3);
    print_matrix_bigint(matrix);

    change_dim_matrix_bigint(&matrix , 4 , 2);
    print_matrix_bigint(matrix);


    matrix_bigint A,B,res_matrix;
    init_matrix_bigint(&A , 3 , 3);
    init_matrix_bigint(&B , 3 , 3);
    init_matrix_bigint(&res_matrix , 3 , 3);

    set_all_coeffs_random_matrix_bigint(A , 10);
    set_all_coeffs_random_matrix_bigint(B , 10);

    add_matrix_bigint(&res_matrix , A , B);

    printf("A : \n");
    print_matrix_bigint(A);

    printf("B : \n");
    print_matrix_bigint(B);

    printf("A+B : \n");
    print_matrix_bigint(res_matrix);

    sub_matrix_bigint(&res_matrix , B , A);

    printf("B-A : \n");
    print_matrix_bigint(res_matrix);

    destroy_matrix_bigint(res_matrix);
    init_matrix_bigint(&res_matrix , 3 , 3);

    mult_matrix_bigint(&res_matrix , A , B);

    printf("A*B (naive) :\n");
    print_matrix_bigint(res_matrix);


    destroy_matrix_bigint(matrix);
    destroy_matrix_bigint(A);
    destroy_matrix_bigint(B);
    destroy_matrix_bigint(res_matrix);

    printf("\n\n*******************************************\n\n");

    matrix_bigQ matrixBigQ , res_matrixBigQ;

    init_matrix_bigQ(&matrixBigQ , 3 , 3);
    init_matrix_bigQ(&res_matrixBigQ , 3 , 3);

    set_coeff_matrix_bigQ_d(matrixBigQ , 0 , 1 , 1 , 2);

    scalar_mult_matrix_bigQ(&matrixBigQ , matrixBigQ , 2);

    change_dim_matrix_bigQ(&matrixBigQ , 2 , 4);

    print_matrix_bigQ(matrixBigQ);

    transpose_matrix_bigQ(&res_matrixBigQ , matrixBigQ);

    print_matrix_bigQ(res_matrixBigQ);

    set_all_coeffs_random_matrix_bigQ(matrixBigQ , 100);


    printf("\n\n*******************************************\n\n");

    matrix_bigQ L , U;
    mpq_t det;
    mpq_init(det);

    matrix_double A_double,L_double , U_double;

    init_matrix_bigQ(&L , 4 , 4);
    init_matrix_bigQ(&U , 4 , 4);

    init_matrix_double(&A_double , 4 , 4);
    init_matrix_double(&L_double , 4 , 4);
    init_matrix_double(&U_double , 4 , 4);


    change_dim_matrix_bigQ(&matrixBigQ , 4 , 4);

    for(unsigned int j=0 ; j<matrixBigQ.nb_col ; j++) {
        set_coeff_matrix_bigQ_d(matrixBigQ , 0 , j , 5 , 1);
        setCoeff_matrix_double(&A_double , 0 , j , 5);
    }

    set_coeff_matrix_bigQ_d(matrixBigQ, 1 , 0 , 10 , 1);
    setCoeff_matrix_double(&A_double , 1 , 0 , 10);

    for(unsigned int j=1 ; j<matrixBigQ.nb_col ; j++) {
        set_coeff_matrix_bigQ_d(matrixBigQ, 1 , j, 16 , 1);
        setCoeff_matrix_double(&A_double , 1 , j , 16);
    }

    set_coeff_matrix_bigQ_d(matrixBigQ, 2 , 0 , 15 , 1);
    set_coeff_matrix_bigQ_d(matrixBigQ, 2 , 1 , 27 , 1);

    setCoeff_matrix_double(&A_double , 2 , 0 , 15);
    setCoeff_matrix_double(&A_double , 2 , 1 , 27);


    for(unsigned int j=2 ; j<matrixBigQ.nb_col ; j++) {
        set_coeff_matrix_bigQ_d(matrixBigQ,2 , j , 34 , 1);
        setCoeff_matrix_double(&A_double , 2 , j , 34);
    }

    set_coeff_matrix_bigQ_d(matrixBigQ, 3 , 0 , 20 , 1);
    set_coeff_matrix_bigQ_d(matrixBigQ, 3 , 1 , 38 , 1);
    set_coeff_matrix_bigQ_d(matrixBigQ, 3 , 2 , 52 , 1);
    set_coeff_matrix_bigQ_d(matrixBigQ, 3 , 3 , 60 , 1);

    setCoeff_matrix_double(&A_double , 3 , 0 , 20);
    setCoeff_matrix_double(&A_double , 3 , 1 , 38);
    setCoeff_matrix_double(&A_double , 3 , 2 , 52);
    setCoeff_matrix_double(&A_double , 3 , 3 , 60);

    LU_decomposition_matrix_bigQ(&L , &U , matrixBigQ);
    LU_decomposition_matrix_double(&L_double , &U_double , A_double);

    printf("A : \n");
    print_matrix_bigQ(matrixBigQ);

    printf("L : \n");
    print_matrix_bigQ(L);

    printf("U : \n");
    print_matrix_bigQ(U);

    print_matrix_double(A_double , "A_double");
    print_matrix_double(L_double , "L_double");
    print_matrix_double(U_double , "U_double");

    determinant_matrix_bigQ(&det , matrixBigQ);
    gmp_printf("\ndet(A) = %Qd\n",det);

    destroy_matrix_bigQ(matrixBigQ);
    destroy_matrix_bigQ(res_matrixBigQ);

    destroy_matrix_double(A_double);
    destroy_matrix_double(L_double);
    destroy_matrix_double(U_double);

    printf("\n\n*******************************************\n\n");

    crt_bigint CRT1 , CRT2, resCRT, resCRT_DAC;

    init_CRT_bigint_d(&CRT1 , 3 , 5);
    init_CRT_bigint_d(&CRT2 , 4 , 7);
    init_CRT_bigint_d(&resCRT , 0 , 0);
    init_CRT_bigint_d(&resCRT_DAC , 0 , 0);

    crt_bigint arg_crt[2];
    arg_crt[0] = CRT1;
    arg_crt[1] = CRT2;


    solve_CRT_bigint(&resCRT , arg_crt ,2);

    print_CRT_bigint_arg(arg_crt , 2);
    printf("\t ***** CRT result *****\n");
    print_CRT_bigint(resCRT);
    printf("\n");

    crt_bigint arg_crt_dac[4];
    for(size_t i=0 ; i<4 ; i++) init_CRT_bigint_d(&arg_crt_dac[i] , 0 , 0);

    mpz_set_d(arg_crt_dac[0].Ai , 3);
    mpz_set_d(arg_crt_dac[0].Ni , 5);

    mpz_set_d(arg_crt_dac[1].Ai , 4);
    mpz_set_d(arg_crt_dac[1].Ni , 7);

    mpz_set_d(arg_crt_dac[2].Ai , 6);
    mpz_set_d(arg_crt_dac[2].Ni , 9);

    mpz_set_d(arg_crt_dac[3].Ai , 11);
    mpz_set_d(arg_crt_dac[3].Ni , 16);

    solve_CRT_bigint_DAC(&resCRT_DAC , arg_crt_dac , 4);

    print_CRT_bigint_arg(arg_crt_dac , 4);
    printf("\t ***** CRT result *****\n");
    print_CRT_bigint(resCRT_DAC);
    printf("\n");

    destroy(CRT1);
    destroy(CRT2);
    destroy(resCRT);
    destroy(resCRT_DAC);

    printf("\n\n*******************************************\n\n");

    pol_bigQ A_bigQ,B_bigQ,Q_bigQ,R_bigQ;

    init_pol_bigQ(&A_bigQ , 5);
    init_pol_bigQ(&B_bigQ , 2);
    init_pol_bigQ(&Q_bigQ , 0);
    init_pol_bigQ(&R_bigQ , 0);

    set_coeff_pol_bigQ_d(A_bigQ , 1 , -2);
    set_coeff_pol_bigQ_d(A_bigQ , 2 , 3);
    set_coeff_pol_bigQ_d(A_bigQ , 3 , -1);
    set_coeff_pol_bigQ_d(A_bigQ , 4 , -1);
    set_coeff_pol_bigQ_d(A_bigQ , 5 , 1);

    set_coeff_pol_bigQ_d(B_bigQ , 0 , 1);
    set_coeff_pol_bigQ_d(B_bigQ , 1 , -1);
    set_coeff_pol_bigQ_d(B_bigQ , 2 , 1);

    euclideDiv_pol_bigQ(&Q_bigQ , &R_bigQ , A_bigQ , B_bigQ);

    printf("\t\t***** Euclidean division bigQ *****\n\n");

    print_pol_bigQ(A_bigQ , "A");
    print_pol_bigQ(B_bigQ , "B");
    print_pol_bigQ(Q_bigQ , "Q");
    print_pol_bigQ(R_bigQ , "R");

    printf("\n\n*******************************************\n\n");

    mpq_t lamdda;
    mpq_init(lamdda);
    mpq_set_si(lamdda , 3 , 4);

    mpq_mult_polbigQ_si(&R_bigQ , A_bigQ , 1 , 2);
    print_pol_bigQ(R_bigQ , "1/2 * A");

    mpq_mult_polbigQ(&R_bigQ , A_bigQ , lamdda);
    print_pol_bigQ(R_bigQ , "3/4 * A");



    printf("\n\n*******************************************\n\n");

    pol_bigQ_value constraints[3];

    init_polbigQ_value_si(&constraints[0] , 0 , 1 , 1 , 1);
    init_polbigQ_value_si(&constraints[1] , 1 , 1 , 4 , 1);
    init_polbigQ_value_si(&constraints[2] , -1 , 1 , 0 , 1);

    print_lagrange_constaints(constraints , 3);

    pol_bigQ M;
    init_pol_bigQ(&M , 2);
    lagrange_interpolation(&M , constraints , 3);
    print_pol_bigQ(M , "Lagrange result");

    printf("\n\n*******************************************\n\n");

    pol_bigQ_value constraints2[4];

    init_polbigQ_value_si(&constraints2[0] , -2 , 1 , 17 , 1);
    init_polbigQ_value_si(&constraints2[1] , -1 , 1 , 13 , 1);
    init_polbigQ_value_si(&constraints2[2] , 1 , 1 , 5 , 1);
    init_polbigQ_value_si(&constraints2[3] , 2 , 1 , 25 , 1);

    print_lagrange_constaints(constraints2 , 4);

    set_all_coeffs_to_pol_bigQ_si(M , 0 , 1);
    lagrange_interpolation(&M , constraints2 , 4);
    print_pol_bigQ(M , "Lagrange results");

    destroy_polbigQ_value(constraints[0]);
    destroy_polbigQ_value(constraints[1]);
    destroy_polbigQ_value(constraints[2]);

    destroy_polbigQ_value(constraints2[0]);
    destroy_polbigQ_value(constraints2[1]);
    destroy_polbigQ_value(constraints2[2]);
    destroy_polbigQ_value(constraints2[3]);

    destroy_pol_bigQ(A_bigQ);
    destroy_pol_bigQ(B_bigQ);
    destroy_pol_bigQ(Q_bigQ);
    destroy_pol_bigQ(R_bigQ);

    printf("\n\n*******************************************\n\n");

    pol_bigQ F_bigQ, G_bigQ;

    init_pol_bigQ(&F_bigQ , 3);
    init_pol_bigQ(&G_bigQ , 2);

    set_coeff_pol_bigQ_d(F_bigQ , 0 , 1);
    set_coeff_pol_bigQ_d(F_bigQ , 1 , 7);
    set_coeff_pol_bigQ_d(F_bigQ , 2 , -2);
    set_coeff_pol_bigQ_d(F_bigQ , 3 , 3);

    set_coeff_pol_bigQ_d(G_bigQ , 0 , 1);
    set_coeff_pol_bigQ_d(G_bigQ , 1 , -1);
    set_coeff_pol_bigQ_d(G_bigQ , 2 , 5);

    print_pol_bigQ(F_bigQ , "F");
    print_pol_bigQ(G_bigQ , "G");

    matrix_bigQ sylvesterFG;
    init_matrix_bigQ(&sylvesterFG , 5 , 5);

    mpq_t resultant;
    mpq_init(resultant);

    sylvester_matrix_bigQ(&sylvesterFG , F_bigQ , G_bigQ);

    printf("Sylvester(F,G) : \n");
    print_matrix_bigQ(sylvesterFG);

    resultant_pol_bigQ(&resultant , F_bigQ , G_bigQ);
    gmp_printf("Res(F,G) : %Qd\n",resultant);

    mpq_set_d(resultant , -1);
    resultant_pol_bigQ_euclidean(&resultant , F_bigQ , G_bigQ);
    gmp_printf("Res(F,G) [Euclidean]: %Qd\n\n",resultant);

    destroy_matrix_bigQ(sylvesterFG);
    destroy_pol_bigQ(F_bigQ);
    destroy_pol_bigQ(G_bigQ);



    mpz_clears(e,f,PM,resM,P,(mpz_t *)NULL);
    mpq_clears(det,lamdda,resultant,(mpq_t*)NULL);

    return 0;
}
