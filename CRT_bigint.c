//
// Created by root on 23/11/18.
//

#include <stdio.h>
#include "CRT_bigint.h"


void init_CRT_bigint(crt_bigint *crt , mpz_t A , mpz_t N) {
    mpz_init_set(crt->Ai , A);
    mpz_init_set(crt->Ni , N);
}

/* ************************************************************************************************************* */

void init_CRT_bigint_d(crt_bigint *crt , int A , int N) {
    mpz_inits(crt->Ai , crt->Ni, (mpz_t*)NULL);
    mpz_set_d(crt->Ai , A);
    mpz_set_d(crt->Ni , N);
//    mpz_init_set_d(crt.Ai , A);
//    mpz_init_set_d(crt.Ni , N);
}

/* ************************************************************************************************************* */

void destroy(crt_bigint crt) {
    mpz_clears(crt.Ni , crt.Ai , (mpz_t*)NULL);
}

/* ************************************************************************************************************* */

void print_CRT_bigint(crt_bigint CRT) {
    gmp_printf("Ai = %Zd\tNi = %Zd\n",CRT.Ai , CRT.Ni);
}

/* ************************************************************************************************************* */

void print_CRT_bigint_arg(crt_bigint *arg , unsigned int nb_arg) {
    printf("\t***** CRT constraints *****\n");
    for(size_t i=0 ; i<nb_arg ; i++) {
       print_CRT_bigint(arg[i]);
    }
}

/* ************************************************************************************************************* */

void solve_CRT_bigint(crt_bigint *res , crt_bigint *arg , unsigned int nb_arg) {

    mpz_t temp, N_on_Ni;
    mpz_inits(temp,N_on_Ni,(mpz_t*)NULL);

    mpz_set_d(res->Ni,1);
    mpz_set_d(res->Ai,0);

    for(size_t i = 0 ; i<nb_arg ; i++) mpz_mul(res->Ni,res->Ni,arg[i].Ni);

    for(size_t i=0 ; i<nb_arg ; i++) {

        mpz_div(N_on_Ni,res->Ni,arg[i].Ni);

        inv_bigint_mod_P(temp,N_on_Ni,arg[i].Ni);

        mpz_mul(temp,temp,N_on_Ni);
        mpz_mul(temp,temp,arg[i].Ai);

        mpz_add(res->Ai,res->Ai,temp);
    }

    mpz_mod(res->Ai,res->Ai,res->Ni);

    mpz_clears(temp,N_on_Ni,(mpz_t*)NULL);
}

/* ************************************************************************************************************* */

void solve_CRT_bigint_DAC(crt_bigint *res , crt_bigint *arg , unsigned int nb_arg) {
    if(log_base_2(nb_arg) - floor(log_base_2(nb_arg)) != 0) {
        perror("nb_arg must be a power of 2");
        return;
    }

    if(nb_arg == 2) {
        solve_CRT_bigint(res , arg , nb_arg);
        return;
    }

    crt_bigint res1 , res2;
    init_CRT_bigint_d(&res1 , 0 , 0);
    init_CRT_bigint_d(&res2 , 0 , 0);

    solve_CRT_bigint_DAC(&res1 , arg , nb_arg/2);
    solve_CRT_bigint_DAC(&res2 , arg+(nb_arg/2) , nb_arg/2);

    crt_bigint crt_arg[2];
    crt_arg[0] = res1;
    crt_arg[1] = res2;

    solve_CRT_bigint(res , crt_arg , 2);

}