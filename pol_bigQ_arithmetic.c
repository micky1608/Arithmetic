//
// Created by root on 30/11/18.
//

#include <stdlib.h>
#include <stdio.h>
#include "pol_bigQ_arithmetic.h"
#include "Util.h"



/* ********************************************************************************************************************** */

/**
 * Allocate the memory for the coefficients
 * @param degre
 */
void init_pol_bigQ(pol_bigQ *polynomial, unsigned int degree_p) {
    if(degree_p >=0) {
        polynomial->degree = degree_p;
        polynomial->coeffs = calloc(degree_p + 1 , sizeof(mpq_t));
        for(int i=0 ; i<= polynomial->degree ; i++) {
            mpq_init(polynomial->coeffs[i]);
            mpq_set_d(polynomial->coeffs[i] , 0);
        }
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
void destroy_pol_bigQ(pol_bigQ polynomial) {
    for(int i=0 ; i<= polynomial.degree ; i++)
        mpq_clear(polynomial.coeffs[i]);

    free(polynomial.coeffs);
}

/* ********************************************************************************************************************** */

/**
 * Change the coefficients of a polynomial from an array containing [ f0, f1 .... ]
 * @param polynomial
 * @param coeffs
 */
void set_all_coeffs_pol_bigQ(pol_bigQ polynomial, mpq_t *coeffs, size_t size_coeffs_array) {
    for(int i=0 ; i<size_coeffs_array ; i++) {
        if (i <= polynomial.degree)
            mpq_set(polynomial.coeffs[i], coeffs[i]);
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
void set_all_coeffs_random_pol_bigQ(pol_bigQ polynomial, unsigned int max) {

    for(int i=0 ; i<=polynomial.degree ; i++) {
        int numerator = (rand() % (2*max)) - max;
        unsigned int denominator;
        do {
            denominator = rand() % max;
        } while(!denominator);
        mpq_set_si(polynomial.coeffs[i],numerator,denominator);
    }
}

/* ********************************************************************************************************************** */

/**
 * Set all the coefficients to the same value
 * @param polynomial
 * @param value
 */
void set_all_coeffs_to_pol_bigQ(pol_bigQ polynomial, mpq_t value) {
    for(int i=0 ; i<=polynomial.degree ; i++)
        mpq_set(polynomial.coeffs[i] , value);
}

/* ********************************************************************************************************************** */

void set_all_coeffs_to_pol_bigQ_si(pol_bigQ polynomial, int numerator , unsigned int denominator) {
    for(int i=0 ; i<=polynomial.degree ; i++)
        mpq_set_si(polynomial.coeffs[i] , numerator , denominator);
}

/* ********************************************************************************************************************** */

/**
 * Change one coefficient in a polynomial
 * @param polynome
 * @param degre_coeff
 * @param newValue
 */
void set_coeff_pol_bigQ(pol_bigQ polynomial, unsigned int degree_coeff, mpq_t newValue) {
    if(degree_coeff <= polynomial.degree)
        mpq_set(polynomial.coeffs[degree_coeff] , newValue);
}

/* ********************************************************************************************************************** */

/**
 * Change one coefficient in a polynomial
 * @param polynome
 * @param degre_coeff
 * @param newValue
 */
void set_coeff_pol_bigQ_d(pol_bigQ polynomial, unsigned int degree_coeff, int newValue) {
    if(degree_coeff <= polynomial.degree) {
        mpq_set_d(polynomial.coeffs[degree_coeff] , newValue);
    }
}

/* ********************************************************************************************************************** */


/**
 * Change the degree of a polynomial
 * If the new degree is greater than the old, the function complete with 0
 * If the new degree is smaller than the old, some coefficients are lost
 * @param polynomial
 */
void change_degre_pol_bigQ(pol_bigQ *polynomial, unsigned int new_degree) {
    if(new_degree == polynomial->degree || new_degree < 0) return;

    mpq_t *new_coeffs = calloc(new_degree+1 , sizeof(mpq_t));

    if(new_degree < polynomial->degree) {
        for(int i=0 ; i<=polynomial->degree ; i++) {
            if(i<=new_degree) {
                mpq_init(new_coeffs[i]);
                mpq_set(new_coeffs[i] , polynomial->coeffs[i]);
            }


            mpq_clear(polynomial->coeffs[i]);
        }
    }
    else {
        for(int i=0 ; i<=new_degree ; i++) {
            if(i<=polynomial->degree) {
                mpq_init(new_coeffs[i]);
                mpq_set(new_coeffs[i] , polynomial->coeffs[i]);
                mpq_clear(polynomial->coeffs[i]);
            }
            else {
                mpq_init(new_coeffs[i]);
                mpq_set_d(new_coeffs[i],0);
            }
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
void copy_pol_bigQ(pol_bigQ *res, pol_bigQ polynomial) {
    change_degre_pol_bigQ(res, polynomial.degree);
    for(int i=0 ; i<=res->degree ; i++)
        mpq_set(res->coeffs[i] , polynomial.coeffs[i]);
}

/* ********************************************************************************************************************** */

void print_pol_bigQ(pol_bigQ polynomial, char *polbigQ_name) {
    mpq_t temp;
    mpq_init(temp);

    printf("%s : ",polbigQ_name);

    if(is_null_polbigQ(polynomial) == TRUE) {
        printf("0\n");
        return;
    }

    if (mpq_cmp_si(polynomial.coeffs[0],0,1) != 0 ) gmp_printf("%Qd",polynomial.coeffs[0]);
    for(int i=1 ; i <= polynomial.degree ; i++) {

        if (mpq_cmp_si(polynomial.coeffs[i],0,1) < 0 )
            printf(" - ");
        else if (mpq_cmp_si(polynomial.coeffs[i],0,1) > 0 )
            printf(" + ");
        else
            continue;

        mpq_abs(temp , polynomial.coeffs[i]);
        printf("(");
        gmp_printf("%Qd",temp);
        printf(" * X^%d)" , i);

    }
    printf("\n");
    mpq_clear(temp);

}

/* ********************************************************************************************************************** */

/**
 * Add two polynomials
 * @param res
 * @param A
 * @param B
 */
void add_pol_bigQ(pol_bigQ *res, pol_bigQ A, pol_bigQ B) {

    if(res->degree != max(A.degree , B.degree))
        change_degre_pol_bigQ(res, max(A.degree, B.degree));

    for(int i=0 ; i<=res->degree ; i++) {
        if(i <= A.degree) {
            if(i <= B.degree) {
                mpq_add(res->coeffs[i] , A.coeffs[i] , B.coeffs[i]);
            }
            else
                set_coeff_pol_bigQ(*res, i, A.coeffs[i]);
        }
        else
            set_coeff_pol_bigQ(*res, i, B.coeffs[i]);
    }
}

/* ********************************************************************************************************************** */

/**
 * Subtract two polynomials
 * @param res
 * @param A
 * @param B
 */
void sub_pol_bigQ(pol_bigQ *res, pol_bigQ A, pol_bigQ B) {

    if(res->degree != max(A.degree , B.degree))
        change_degre_pol_bigQ(res, max(A.degree, B.degree));

    for(int i=0 ; i<=res->degree ; i++) {
        if(i <= A.degree) {
            if(i <= B.degree) {
                mpq_sub(res->coeffs[i] , A.coeffs[i] , B.coeffs[i]);
            }
            else
                set_coeff_pol_bigQ(*res, i, A.coeffs[i]);
        }
        else
            mpq_neg(*(res->coeffs+i) , B.coeffs[i]);
    }

}

/* ********************************************************************************************************************** */

/**
 * Naive multiplication of two polynomials
 * @param res
 * @param A
 * @param B
 */
void mult_pol_bigQ(pol_bigQ *res, pol_bigQ A, pol_bigQ B) {
    change_degre_pol_bigQ(res, A.degree + B.degree);
    mpq_t temp;
    mpq_init(temp);

    for(int i=0 ; i<=res->degree ; i++) {
        mpq_set_d(res->coeffs[i] , 0);

        int index_A = 0;
        int index_B = i;

        while (index_A <= i && index_B >=0) {

            if(index_A <= A.degree && index_B <= B.degree) {
                mpq_mul(temp , A.coeffs[index_A] , B.coeffs[index_B]);
                mpq_add(res->coeffs[i] , res->coeffs[i] , temp);
            }
            index_A++;
            index_B--;
        }
    }

    mpq_clear(temp);
}

/* ********************************************************************************************************************** */

void mpq_mult_polbigQ(pol_bigQ *res, pol_bigQ A, mpq_t lambda) {
    change_degre_pol_bigQ(res , A.degree);

    for(int i=0 ; i<=res->degree ; i++) {
            mpq_mul(res->coeffs[i] , A.coeffs[i] , lambda);
    }
}

/* ********************************************************************************************************************** */

void mpq_mult_polbigQ_si(pol_bigQ *res, pol_bigQ A, int lambda_num, unsigned int lambda_den) {
    change_degre_pol_bigQ(res , A.degree);

    mpz_t numerator, denominator,gcd;
    mpz_inits(numerator , denominator , gcd , NULL);

    for(int i=0 ; i<=res->degree ; i++) {
        mpq_get_num(numerator , A.coeffs[i]);
        mpq_get_den(denominator , A.coeffs[i]);

        mpz_mul_si(numerator , numerator , lambda_num);
        mpz_mul_si(denominator , denominator , lambda_den);

        mpz_gcd(gcd , numerator , denominator);
        mpz_div(numerator , numerator , gcd);
        mpz_div(denominator , denominator , gcd);

        mpq_set_num(res->coeffs[i] , numerator);
        mpq_set_den(res->coeffs[i] , denominator);

    }

    mpz_clears(numerator , denominator , gcd , NULL);
}

/* ********************************************************************************************************************** */

/**
 * Degree A = Degree B = 2^p -1
 * @param res is a polynomial of degree 2*degree(A) ! THE MEMORY MUST BE ALLOCATED BEFORE CALLING THIS FUNCTION
 * @param A
 * @param B
 */
void karatsuba_pol_bigQ(pol_bigQ *res, pol_bigQ A, pol_bigQ B) {


    if(A.degree != B.degree) {
        perror("Karatsuba with different degrees");
        return;
    }

    if(log_base_2(A.degree + 1) - floor(log_base_2(A.degree + 1)) != 0) {
        perror("Karatsuba degree != 2^p - 1");
        return;
    }


    if(A.degree == 1) {
        mpq_mul(res->coeffs[0] , A.coeffs[0] , B.coeffs[0]);
        mpq_mul(res->coeffs[2] , A.coeffs[1] , B.coeffs[1]);

        mpq_t sum_coeffs_A , sum_coeffs_B;

        mpq_inits(sum_coeffs_A , sum_coeffs_B,(mpq_t *)NULL);

        mpq_add(sum_coeffs_A , A.coeffs[0] , A.coeffs[1]);
        mpq_add(sum_coeffs_B , B.coeffs[0] , B.coeffs[1]);

        mpq_mul(res->coeffs[1] , sum_coeffs_A , sum_coeffs_B);
        mpq_sub(res->coeffs[1] , res->coeffs[1] , res->coeffs[0]);
        mpq_sub(res->coeffs[1] , res->coeffs[1] , res->coeffs[2]);

        mpq_clears(sum_coeffs_A , sum_coeffs_B , (mpq_t *)NULL);

        return;
    }

    size_t p = (size_t)(log_base_2(A.degree + 1));

    pol_bigQ A_under, A_upper , B_under , B_upper , sum_A_under_upper , sum_B_under_upper , H0 , H1 , H2;

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

    init_pol_bigQ(&H1, H0.degree);
    init_pol_bigQ(&sum_A_under_upper, A_under.degree);
    init_pol_bigQ(&sum_B_under_upper, B_under.degree);

    add_pol_bigQ(&sum_A_under_upper, A_under, A_upper);
    add_pol_bigQ(&sum_B_under_upper, B_under, B_upper);

    karatsuba_pol_bigQ(&H0, A_under, B_under);
    karatsuba_pol_bigQ(&H2, A_upper, B_upper);
    karatsuba_pol_bigQ(&H1, sum_A_under_upper, sum_B_under_upper);

    sub_pol_bigQ(&H1, H1, H0);
    sub_pol_bigQ(&H1, H1, H2);


    for(int i=0 , index_begin = (int)pow(2 , (p-1)) ; i<= H1.degree ; i++)
        mpq_add(res->coeffs[index_begin + i] , res->coeffs[index_begin + i] , H1.coeffs[i]);

    destroy_pol_bigQ(H1);
    destroy_pol_bigQ(sum_A_under_upper);
    destroy_pol_bigQ(sum_B_under_upper);

}

/* ********************************************************************************************************************** */

void euclideDiv_pol_bigQ(pol_bigQ *Q , pol_bigQ *R , pol_bigQ A , pol_bigQ B) {
    if(A.degree <= B.degree) {
        perror("Euclidean division : A.degree <= B.degree !!");
        return;
    }

    mpq_t a,b, a_on_b, a_on_b_neg;

    mpq_inits(a,b, a_on_b , a_on_b_neg , (mpq_t *)NULL);
    mpq_set(b, B.coeffs[B.degree]);

    pol_bigQ temp , temp2;

    init_pol_bigQ(&temp, A.degree - B.degree);
    init_pol_bigQ(&temp2, A.degree);

    change_degre_pol_bigQ(Q, A.degree - B.degree);
    copy_pol_bigQ(R, A);
    set_all_coeffs_to_pol_bigQ_si(*Q, 0,1);


    while(R->degree >= B.degree) {
        //for(int k=0 ; k<1 ; k++) {
        mpq_set(a,R->coeffs[R->degree]);
        mpq_div(a_on_b , a , b);
        mpq_neg(a_on_b_neg , a_on_b);

        set_all_coeffs_to_pol_bigQ_si(temp, 0,1);
        if(temp.degree != R->degree-B.degree) change_degre_pol_bigQ(&temp, R->degree - B.degree);
        set_coeff_pol_bigQ(temp, R->degree - B.degree, a_on_b);

        add_pol_bigQ(Q, *Q, temp);

        set_coeff_pol_bigQ(temp, R->degree - B.degree, a_on_b_neg);

        if(temp2.degree != R->degree) change_degre_pol_bigQ(&temp2 , R->degree);

        mult_pol_bigQ(&temp2, temp, B);

        add_pol_bigQ(R, *R, temp2);

        // update the degree of R
        unsigned int i;

        for(i=R->degree ; mpq_cmp_si(R->coeffs[i],0,1) == 0 && i>0; --i);

        change_degre_pol_bigQ(R , i);

    }

    mpq_clears(a, b, a_on_b, a_on_b_neg, (mpq_t*)NULL);
    destroy_pol_bigQ(temp);
    destroy_pol_bigQ(temp2);

}

/* ********************************************************************************************************************** */

bool is_null_polbigQ(pol_bigQ A) {
    for(int i=0 ; i<=A.degree ; i++) {
        if(mpq_cmp_si(A.coeffs[i] , 0 , 1) != 0)
            return FALSE;
    }
    return TRUE;
}

/* ********************************************************************************************************************** */

void sylvester_matrix_bigQ(matrix_bigQ *sylvester , pol_bigQ F , pol_bigQ G) {
    unsigned int p = F.degree , q=G.degree;
    change_dim_matrix_bigQ(sylvester , p+q , p+q);

    int j;
    for(int i=0 ; i<q ; i++) {
        for(j=0 ; j<i ; j++) mpq_set_d(sylvester->values[i*sylvester->nb_col+j] , 0);
        for(int k=p ; k>=0 ; k--) {
            mpq_set(sylvester->values[i*sylvester->nb_col+j] , F.coeffs[k]);
            j++;
        }
        while(j<p+q) {
            mpq_set_d(sylvester->values[i*sylvester->nb_col+j] , 0);
            j++;
        }
    }


    for(int i=q ; i<p+q ; i++) {
        for(j=0 ; j<i-q ; j++) mpq_set_d(sylvester->values[i*sylvester->nb_col+j] , 0);
        for(int k=q ; k>=0 ; k--) {
            mpq_set(sylvester->values[i*sylvester->nb_col+j] , G.coeffs[k]);
            j++;
        }
        while(j<p+q) {
            mpq_set_d(sylvester->values[i*sylvester->nb_col+j] , 0);
            j++;
        }
    }
}

/* ********************************************************************************************************************** */

void resultant_pol_bigQ(mpq_t *resultant , pol_bigQ F , pol_bigQ G) {
    matrix_bigQ sylvester;
    init_matrix_bigQ(&sylvester , F.degree+G.degree , F.degree+G.degree);
    sylvester_matrix_bigQ(&sylvester , F , G);
    determinant_matrix_bigQ(resultant , sylvester);
    destroy_matrix_bigQ(sylvester);
}


