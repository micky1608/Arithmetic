//
// Created by root on 12/10/18.
//

#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include "polynomial_arithmetic.h"


/**
 * Allocate the memory for the coefficients
 * @param degre
 */
void init_Pol_M(Pol_M *polynomial , unsigned int degree_p) {
    polynomial->degree = degree_p;
    polynomial->coeffs = calloc(degree_p + 1 , sizeof(mpz_t));
    for(int i=0 ; i<= polynomial->degree ; i++) mpz_init(polynomial->coeffs[i]);
}

/**
 * Free the memery allocated for a polynomial
 * @param polynomial
 */
void destroy_Pol_M(Pol_M polynomial) {
    for(int i=0 ; i<= polynomial.degree ; i++) mpz_clear(polynomial.coeffs[i]);
    free(polynomial.coeffs);
}

/**
 * Change all the coefficients of a polynomial from an array containing [ f0, f1 .... ]
 * @param polynomial
 * @param coeffs
 */
void set_all_coeffs_Pol_M(Pol_M polynomial , mpz_t *coeffs) {
    //TODO
}
/**
 * Define all the coefficients of a polynomial between -max and max
 * @param polynomial
 * @param max
 */
void set_all_coeffs_random_Pol_M(Pol_M polynomial , unsigned int max) {
    srand(time(NULL));

    for(int i=0 ; i<=polynomial.degree ; i++) {
        int fi_i = (rand() % (2*max)) - max;
        mpz_t fi;
        mpz_set_d(polynomial.coeffs[i],fi_i);
    }

}

/**
 * Change one coefficient in a polynomial
 * @param polynome
 * @param degre_coeff
 * @param newValue
 */
void set_coeff_Pol_M(Pol_M polynomial , unsigned int degree_coeff , mpz_t newValue) {
    if(degree_coeff <= polynomial.degree)
        mpz_set(polynomial.coeffs[degree_coeff] , newValue);
}


/**
 * Change the degree of a polynomial
 * If the new degree is greater than the old, the function complete with 0
 * If the new degree is smaller than the old, some coefficients are lost
 * @param polynomial
 */
void change_degre_Pol_M(Pol_M *polynomial, unsigned int new_degree) {
    if(new_degree == polynomial->degree) return;

    if(new_degree < polynomial->degree) {
        for(int i=new_degree + 1  ; i<=polynomial->degree ; i++) mpz_clear(polynomial->coeffs[i]);
        polynomial->coeffs = realloc(polynomial->coeffs , new_degree);
    }
    else {
        polynomial->coeffs = realloc(polynomial->coeffs , new_degree + 1);
        for(int i=polynomial->degree + 1 ; i <= new_degree ; i++) {
            mpz_init(polynomial->coeffs[i]);
            mpz_set_ui(polynomial->coeffs[i] , 0);
        }
    }

    polynomial->degree = new_degree;
}

void print_Pol_m(Pol_M polynomial) {
    mpz_out_str(stdout,10,polynomial.coeffs[0]);
    printf(" + ");
    for(int i=1 ; i <= polynomial.degree ; i++) {
        mpz_out_str(stdout , 10 , polynomial.coeffs[i]);
        printf(" * X^%d" , i);
        if(i!=polynomial.degree)
            printf(" + ");
    }
    printf("\n");

}

