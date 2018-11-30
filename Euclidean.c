//
// Created by root on 29/10/18.
//

#include "Euclidean.h"



/* ********************************************************************************************************************** */

/**
 * Euclidean division of A by B with A and B two polynomials
 * Q and R must have mpf_t coefficients
 * @param Q
 * @param R
 * @param A
 * @param B
 */
void euclideDiv_pol_bignumber(pol_bigfloat *Q , pol_bigfloat *R , pol_bigint A , pol_bigint B) {

    if(A.degree <= B.degree) {
        perror("Euclidean division : A.degree <= B.degree !!");
        return;
    }

    mpf_t a,b, a_on_b, a_on_b_neg, zero;

    mpf_inits(a,b, a_on_b , a_on_b_neg , zero, (mpf_t *)NULL);
    mpf_set_z(b, B.coeffs[B.degree]);
    mpf_set_d(zero , 0);

    pol_bigfloat temp , temp2 , A_float, B_float;

    init_pol_bigfloat(&temp, A.degree - B.degree);
    init_pol_bigfloat(&temp2, A.degree);
    init_pol_bigfloat(&A_float , 0);
    init_pol_bigfloat(&B_float , 0);

    pol_bigint_to_pol_bigfloat(&A_float , A);
    pol_bigint_to_pol_bigfloat(&B_float , B);

    change_degre_pol_bigfloat(Q, A.degree - B.degree);
    copy_pol_bigfloat(R, A_float);
    set_all_coeffs_to_pol_bigfloat(*Q, zero);


    while(R->degree >= B.degree) {
    //for(int k=0 ; k<1 ; k++) {
        mpf_set(a,R->coeffs[R->degree]);
        mpf_div(a_on_b , a , b);
        mpf_neg(a_on_b_neg , a_on_b);

        set_all_coeffs_to_pol_bigfloat(temp, zero);
        if(temp.degree != R->degree-B.degree) change_degre_pol_bigfloat(&temp, R->degree - B.degree);
        set_coeff_pol_bigfloat(temp, R->degree - B.degree, a_on_b);

        add_pol_bigfloat(Q, *Q, temp);

        set_coeff_pol_bigfloat(temp, R->degree - B.degree, a_on_b_neg);

        if(temp2.degree != R->degree) change_degre_pol_bigfloat(&temp2 , R->degree);

        mult_pol_bigfloat(&temp2, temp, B_float);

        add_pol_bigfloat(R, *R, temp2);

        // update the degree of R
        unsigned int i;

        for(i=R->degree ; mpf_cmp_d(R->coeffs[i],0) == 0 ; --i);

        change_degre_pol_bigfloat(R , i);

    }


    mpf_clears(a, b, a_on_b, a_on_b_neg, zero, (mpf_t*)NULL);
    destroy_pol_bigfloat(temp);
    destroy_pol_bigfloat(temp2);
    destroy_pol_bigfloat(A_float);
    destroy_pol_bigfloat(B_float);

}
