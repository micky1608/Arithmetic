//
// Created by root on 12/10/18.
//


#include "pol_bigint_arithmetic.h"
#include "arithmetic.h"


/* ********************************************************************************************************************** */

/**
 * Allocate the memory for the coefficients
 * @param degre
 */
void init_pol_bigint(pol_bigint *polynomial, unsigned int degree_p) {
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

/* ********************************************************************************************************************** */

/**
 * Free the memery allocated for a polynomial
 * @param polynomial
 */
void destroy_pol_bigint(pol_bigint polynomial) {
    for(int i=0 ; i<= polynomial.degree ; i++)
        mpz_clear(polynomial.coeffs[i]);
    free(polynomial.coeffs);
}

/* ********************************************************************************************************************** */

/**
 * Change the coefficients of a polynomial from an array containing [ f0, f1 .... ]
 * @param polynomial
 * @param coeffs
 */
void set_all_coeffs_pol_bigint(pol_bigint polynomial, mpz_t *coeffs, size_t size_coeffs_array) {
    for(int i=0 ; i<size_coeffs_array ; i++) {
        if (i <= polynomial.degree)
            mpz_set(polynomial.coeffs[i], coeffs[i]);
        else
            break;
    }
}

/* ********************************************************************************************************************** */

/**
 * Define all the coefficients of a polynomial between -max and max
 * @param polynomial
 * @param max
 */
void set_all_coeffs_random_pol_bigint(pol_bigint polynomial, unsigned int max) {

    for(int i=0 ; i<=polynomial.degree ; i++) {
        int fi_i = (rand() % (2*max)) - max;
        while(i==polynomial.degree && fi_i==0) {
            fi_i = (rand() % (2*max)) - max;
        }
        mpz_set_d(polynomial.coeffs[i],fi_i);
    }
}

/* ********************************************************************************************************************** */

/**
 * Set all the coefficients to the same value
 * @param polynomial
 * @param value
 */
void set_all_coeffs_to_pol_bigint(pol_bigint polynomial, mpz_t value) {
    for(int i=0 ; i<=polynomial.degree ; i++)
        mpz_set(polynomial.coeffs[i] , value);

}

/* ********************************************************************************************************************** */

/**
 * Change one coefficient in a polynomial
 * @param polynome
 * @param degre_coeff
 * @param newValue
 */
void set_coeff_pol_bigint(pol_bigint polynomial, unsigned int degree_coeff, mpz_t newValue) {
    if(degree_coeff <= polynomial.degree)
        mpz_set(polynomial.coeffs[degree_coeff] , newValue);
}

/* ********************************************************************************************************************** */

/**
 * Change one coefficient in a polynomial
 * @param polynome
 * @param degre_coeff
 * @param newValue
 */
void set_coeff_pol_bigint_d(pol_bigint polynomial, unsigned int degree_coeff, int newValue) {
    if(degree_coeff <= polynomial.degree) {
        mpz_set_d(polynomial.coeffs[degree_coeff] , newValue);
    }

}

/* ********************************************************************************************************************** */


/**
 * Change the degree of a polynomial
 * If the new degree is greater than the old, the function complete with 0
 * If the new degree is smaller than the old, some coefficients are lost
 * @param polynomial
 */
void change_degre_pol_bigint(pol_bigint *polynomial, unsigned int new_degree) {
    if(new_degree == polynomial->degree || new_degree < 0) return;

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

/* ********************************************************************************************************************** */

/**
 * Copy the polynomial src in res
 * @param res
 * @param polynomial
 */
void copy_pol_bigint(pol_bigint *res, pol_bigint polynomial) {
    change_degre_pol_bigint(res, polynomial.degree);
    for(int i=0 ; i<=res->degree ; i++)
        mpz_set(res->coeffs[i] , polynomial.coeffs[i]);
}

/* ********************************************************************************************************************** */

void print_pol_bigint(pol_bigint polynomial) {
    mpz_t temp;
    mpz_init(temp);


    if (mpz_cmp_ui(polynomial.coeffs[0],0) != 0 ) gmp_printf("%Zd",polynomial.coeffs[0]);
    for(int i=1 ; i <= polynomial.degree ; i++) {

        if (mpz_cmp_ui(polynomial.coeffs[i],0) < 0 )
            printf(" - ");
        else if (mpz_cmp_ui(polynomial.coeffs[i],0) > 0 )
            printf(" + ");
        else
            continue;

        mpz_abs(temp , polynomial.coeffs[i]);
        printf("(");
        gmp_printf("%Zd",temp);
        printf(" * X^%d)" , i);

    }
    printf("\n");
    mpz_clear(temp);

}

/* ********************************************************************************************************************** */

/**
 * Add two polynomials
 * @param res
 * @param A
 * @param B
 */
void add_pol_bigint(pol_bigint *res, pol_bigint A, pol_bigint B) {

    if(res->degree != max(A.degree , B.degree))
        change_degre_pol_bigint(res, max(A.degree, B.degree));

    for(int i=0 ; i<=res->degree ; i++) {
        if(i <= A.degree) {
            if(i <= B.degree) {
               mpz_add(res->coeffs[i] , A.coeffs[i] , B.coeffs[i]);
            }
            else
                set_coeff_pol_bigint(*res, i, A.coeffs[i]);
        }
        else
            set_coeff_pol_bigint(*res, i, B.coeffs[i]);
    }
}

/* ********************************************************************************************************************** */

/**
 * Subtract two polynomials
 * @param res
 * @param A
 * @param B
 */
void sub_pol_bigint(pol_bigint *res, pol_bigint A, pol_bigint B) {

    if(res->degree != max(A.degree , B.degree))
        change_degre_pol_bigint(res, max(A.degree, B.degree));

    for(int i=0 ; i<=res->degree ; i++) {
        if(i <= A.degree) {
            if(i <= B.degree) {
                mpz_sub(res->coeffs[i] , A.coeffs[i] , B.coeffs[i]);
            }
            else
                set_coeff_pol_bigint(*res, i, A.coeffs[i]);
        }
        else
            mpz_neg(*(res->coeffs+i) , B.coeffs[i]);
    }

}

/* ********************************************************************************************************************** */

/**
 * Naive multiplication of two polynomials
 * @param res
 * @param A
 * @param B
 */
void mult_pol_bigint(pol_bigint *res, pol_bigint A, pol_bigint B) {
    change_degre_pol_bigint(res, A.degree + B.degree);
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

/* ********************************************************************************************************************** */

/**
 * Degree A = Degree B = 2^p -1
 * @param res is a polynomial of degree 2*degree(A) ! THE MEMORY MUST BE ALLOCATED BEFORE CALLING THIS FUNCTION
 * @param A
 * @param B
 */
void karatsuba_pol_bigint(pol_bigint *res, pol_bigint A, pol_bigint B) {


    if(A.degree != B.degree) {
        perror("Karatsuba with different degrees");
        return;
    }

    if(log_base_2(A.degree + 1) - floor(log_base_2(A.degree + 1)) != 0) {
        perror("Karatsuba degree != 2^p - 1");
        return;
    }


    if(A.degree == 1) {
        mpz_mul(res->coeffs[0] , A.coeffs[0] , B.coeffs[0]);
        mpz_mul(res->coeffs[2] , A.coeffs[1] , B.coeffs[1]);

        mpz_t sum_coeffs_A , sum_coeffs_B;

        mpz_inits(sum_coeffs_A , sum_coeffs_B,(mpz_t *)NULL);

        mpz_add(sum_coeffs_A , A.coeffs[0] , A.coeffs[1]);
        mpz_add(sum_coeffs_B , B.coeffs[0] , B.coeffs[1]);

        mpz_mul(res->coeffs[1] , sum_coeffs_A , sum_coeffs_B);
        mpz_sub(res->coeffs[1] , res->coeffs[1] , res->coeffs[0]);
        mpz_sub(res->coeffs[1] , res->coeffs[1] , res->coeffs[2]);

        mpz_clears(sum_coeffs_A , sum_coeffs_B , (mpz_t *)NULL);

        return;
    }

    size_t p = (size_t)(log_base_2(A.degree + 1));

    pol_bigint A_under, A_upper , B_under , B_upper , sum_A_under_upper , sum_B_under_upper , H0 , H1 , H2;

    A_under.degree = (unsigned int)(pow((double)2,(double)(p-1))) - 1; // 2^(p-1) - 1
    A_upper.degree = A_under.degree;
    B_under.degree = A_under.degree;
    B_upper.degree = A_under.degree;
    H0.degree = 2*A_under.degree;
    H2.degree = H0.degree;

    A_under.coeffs = A.coeffs;
    A_upper.coeffs = &A.coeffs[A_under.degree + 1];
    B_under.coeffs = B.coeffs;
    B_upper.coeffs = &B.coeffs[B_under.degree + 1];

    H0.coeffs = res->coeffs;
    H2.coeffs = &res->coeffs[res->degree - H2.degree];

    init_pol_bigint(&H1, H0.degree);
    init_pol_bigint(&sum_A_under_upper, A_under.degree);
    init_pol_bigint(&sum_B_under_upper, B_under.degree);

    add_pol_bigint(&sum_A_under_upper, A_under, A_upper);
    add_pol_bigint(&sum_B_under_upper, B_under, B_upper);

    karatsuba_pol_bigint(&H0, A_under, B_under);
    karatsuba_pol_bigint(&H2, A_upper, B_upper);
    karatsuba_pol_bigint(&H1, sum_A_under_upper, sum_B_under_upper);

    sub_pol_bigint(&H1, H1, H0);
    sub_pol_bigint(&H1, H1, H2);


    for(int i=0 , index_begin = (int)pow(2 , (p-1)) ; i<= H1.degree ; i++)
        mpz_add(res->coeffs[index_begin + i] , res->coeffs[index_begin + i] , H1.coeffs[i]);

    destroy_pol_bigint(H1);
    destroy_pol_bigint(sum_A_under_upper);
    destroy_pol_bigint(sum_B_under_upper);

}



