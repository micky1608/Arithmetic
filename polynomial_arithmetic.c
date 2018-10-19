//
// Created by root on 12/10/18.
//

#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include "polynomial_arithmetic.h"



/**
 * Allocate the memory for the coefficients
 * @param degre
 */
void init_Pol_M(Pol_M *polynomial , unsigned int degree_p) {
    if(degree_p >=0) {
        polynomial->degree = degree_p;
        polynomial->coeffs = calloc(degree_p + 1 , sizeof(mpz_t));
        for(int i=0 ; i<= polynomial->degree ; i++) mpz_init(polynomial->coeffs[i]);
    }
    else {
        perror("Try to create a polynomial with negative degree !!");
        exit(-1);
    }

}

/**
 * Free the memery allocated for a polynomial
 * @param polynomial
 */
void destroy_Pol_M(Pol_M polynomial) {
    for(int i=0 ; i<= polynomial.degree ; i++)
        mpz_clear(polynomial.coeffs[i]);
    free(polynomial.coeffs);
}

/**
 * Change the coefficients of a polynomial from an array containing [ f0, f1 .... ]
 * @param polynomial
 * @param coeffs
 */
void set_all_coeffs_Pol_M(Pol_M polynomial , mpz_t *coeffs , size_t size_coeffs_array) {
    for(int i=0 ; i<size_coeffs_array ; i++) {
        if (i <= polynomial.degree)
            mpz_set(polynomial.coeffs[i], coeffs[i]);
        else
            break;
    }
}
/**
 * Define all the coefficients of a polynomial between -max and max
 * @param polynomial
 * @param max
 */
void set_all_coeffs_random_Pol_M(Pol_M polynomial , unsigned int max) {

    for(int i=0 ; i<=polynomial.degree ; i++) {
        int fi_i = (rand() % (2*max)) - max;
        while(i==polynomial.degree && fi_i==0) {
            fi_i = (rand() % (2*max)) - max;
        }
        mpz_set_d(polynomial.coeffs[i],fi_i);
    }


}

/**
 * Set all the coefficients to the same value
 * @param polynomial
 * @param value
 */
void set_all_coeffs_to_Pol_M(Pol_M polynomial , mpz_t value) {
    for(int i=0 ; i<=polynomial.degree ; i++)
        mpz_set(polynomial.coeffs[i] , value);

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

    mpz_t *new_coeffs = calloc(new_degree+1 , sizeof(mpz_t));

    if(new_degree < polynomial->degree) {
        for(int i=0 ; i<=polynomial->degree ; i++) {
            if(i<=new_degree)
                mpz_init_set(new_coeffs[i] , polynomial->coeffs[i]);

            mpz_clear(polynomial->coeffs[i]);
        }
    }
    else {
            for(int i=0 ; i<=new_degree ; i++) {
                if(i<=polynomial->degree) {
                    mpz_init_set(new_coeffs[i] , polynomial->coeffs[i]);
                    mpz_clear(polynomial->coeffs[i]);
                }
                else
                    mpz_init_set_d(new_coeffs[i],0);
            }
    }

    free(polynomial->coeffs);
    polynomial->coeffs = new_coeffs;
    polynomial->degree = new_degree;
}

/**
 * Copy the polynomial src in res
 * @param res
 * @param polynomial
 */
void copy_Pol_M(Pol_M *res , Pol_M polynomial) {
    change_degre_Pol_M(res , polynomial.degree);
    for(int i=0 ; i<=res->degree ; i++)
        mpz_set(res->coeffs[i] , polynomial.coeffs[i]);
}

void print_Pol_m(Pol_M polynomial) {
    mpz_t temp;
    mpz_init(temp);


    mpz_out_str(stdout,10,polynomial.coeffs[0]);
    for(int i=1 ; i <= polynomial.degree ; i++) {

        if (mpz_cmp_ui(polynomial.coeffs[i],0) < 0 )
            printf(" - ");
        else
            printf(" + ");

        mpz_abs(temp , polynomial.coeffs[i]);
        printf("(");
        mpz_out_str(stdout , 10 , temp);
        printf(" * X^%d)" , i);

    }
    printf("\n");
    mpz_clear(temp);

}


/**
 * Add two polynomials
 * @param res
 * @param A
 * @param B
 */
void add_Pol_M(Pol_M *res , Pol_M A , Pol_M B) {
    change_degre_Pol_M(res , max(A.degree,B.degree));

    for(int i=0 ; i<=res->degree ; i++) {
        if(i <= A.degree) {
            if(i <= B.degree) {
               mpz_add(res->coeffs[i] , A.coeffs[i] , B.coeffs[i]);
            }
            else
                set_coeff_Pol_M(*res , i , A.coeffs[i]);

        }
        else
            set_coeff_Pol_M(*res , i , B.coeffs[i]);
    }
}

/**
 * Naive multiplication of two polynomials
 * @param res
 * @param A
 * @param B
 */
void mult_Pol_M(Pol_M *res , Pol_M A , Pol_M B) {
    change_degre_Pol_M(res , A.degree+B.degree);
    mpz_t temp;
    mpz_init(temp);

    for(int i=0 ; i<=res->degree ; i++) {
        mpz_set_d(res->coeffs[i] , 0);

        int index_A = 0;
        int index_B = i;

        while (index_A <= i && index_B >=0) {

            if(index_A <= A.degree && index_B <= B.degree) {
                mpz_mul(temp , A.coeffs[index_A] , B.coeffs[index_B]);
                mpz_add(res->coeffs[i] , res->coeffs[i] , temp);
            }
            index_A++;
            index_B--;
        }
    }

    mpz_clear(temp);
}

/**
 * Euclidean division of A by B with A and B two polynomials in mpz/P.mpz
 * @param Q
 * @param R
 * @param A
 * @param B
 */
void euclide_div_Pol_M_Mod_P(Pol_M *Q , Pol_M *R , Pol_M A , Pol_M B , mpz_t P) {
    if(A.degree <= B.degree) {
        perror("Euclidean division : A.degree <= B.degree !!");
        return;
    }

    mpz_t a,b,zero, a_on_b, a_on_b_neg;
    mpz_inits(a, b,zero, a_on_b , a_on_b_neg,0);
    mpz_set(b, B.coeffs[B.degree]);
    mpz_set_d(zero , 0);

    Pol_M temp , temp2;
    init_Pol_M(&temp , A.degree-B.degree);
    init_Pol_M(&temp2 , A.degree-B.degree);

    change_degre_Pol_M(Q , A.degree-B.degree);
    copy_Pol_M(R , A);
    set_all_coeffs_to_Pol_M(*Q , zero);

    while(R->degree > B.degree) {
        mpz_set(a,R->coeffs[R->degree]);
        mpz_fdiv_q(a_on_b , a , b);
        mpz_neg(a_on_b_neg , a_on_b);

        set_all_coeffs_to_Pol_M(temp , zero);
        if(temp.degree != R->degree-B.degree) change_degre_Pol_M(&temp , R->degree-B.degree);
        set_coeff_Pol_M(temp , R->degree-B.degree , a_on_b);

        add_Pol_M(Q , *Q , temp);

        set_coeff_Pol_M(temp , R->degree-B.degree , a_on_b_neg);

        mult_Pol_M(&temp2 , temp , B);

        printf("R : ");
        print_Pol_m(*R);
        printf("\n");

        printf("temp : ");
        print_Pol_m(temp);
        printf("\n");

        printf("temp2 : ");
        print_Pol_m(temp2);
        printf("\n");

        add_Pol_M(R , *R , temp2);
    }

    mpz_clears(a,b,zero,a_on_b,a_on_b_neg,0);


}



