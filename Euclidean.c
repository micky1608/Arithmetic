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

        for(i=R->degree ; mpf_cmp_d(R->coeffs[i],0) == 0 && i>0; --i);

        change_degre_pol_bigfloat(R , i);

    }


    mpf_clears(a, b, a_on_b, a_on_b_neg, zero, (mpf_t*)NULL);
    destroy_pol_bigfloat(temp);
    destroy_pol_bigfloat(temp2);
    destroy_pol_bigfloat(A_float);
    destroy_pol_bigfloat(B_float);
}

/* ********************************************************************************************************************** */

void halfGCD(matrix_pol *Mgcd , pol A , pol B) {
    if(B.degree >= A.degree) {
        perror("HalfGCD B degree must be smaller");
        exit(EXIT_FAILURE);
    }

    int n = A.degree , m = (int)ceil((double)n/2);
    
    if(B.degree < m) {
        identity_matrix_pol(Mgcd , 2);
        return;
    }

    init_matrix_pol(Mgcd , 2 , 2);

    pol x_pow_m , f , g , r , Q , Qneg , x_pow_l , b , c;
    matrix_pol M , AB , ABprime , Mprime , BCprime , M_ABprime , Msecond , temp;

    init_pol(&x_pow_m , (unsigned int)m);
    set_coeff_pol(x_pow_m , (unsigned int)m , 1);

    init_pol(&f , m);
    init_pol(&g , m);
    init_pol(&r , m);
    init_pol(&Q , A.degree - B.degree);
    init_pol(&Qneg , A.degree-B.degree);

    init_matrix_pol(&AB , 2 , 1);
    init_matrix_pol(&ABprime , 2 , 1);
    init_matrix_pol(&Mprime , 2 , 1);
    init_matrix_pol(&BCprime , 2 , 1);
    init_matrix_pol(&M_ABprime , 2 , 2);
    init_matrix_pol(&temp , 2 , 2);

    setCoeff_matrix_pol(&AB , 0 , 0 , A);
    setCoeff_matrix_pol(&AB , 1 , 0 , B);

    euclide_div_pol(&f , &r , A , x_pow_m);
    euclide_div_pol(&g , &r , B , x_pow_m);

    print_pol(f , "f");
    print_pol(g , "g");

    halfGCD(&M , f , g); // recursive call

    mul_matrix_pol(&ABprime , M , AB);

    if(ABprime.values[1].degree < m) {
        copy_matrix_pol(Mgcd , M);
        return;
    }

    euclide_div_pol(&Q , &r , ABprime.values[0] , ABprime.values[1]);
    scalar_mult_pol(&Qneg , Q , -1);
    set_coeff_constant_matrix_pol(&M_ABprime , 0 , 1 , 1);
    set_coeff_constant_matrix_pol(&M_ABprime , 1 , 0 , 1);
    setCoeff_matrix_pol(&M_ABprime , 1 , 1 , Qneg);

    mul_matrix_pol(&BCprime , M_ABprime , ABprime);

    int l = 2*m - ABprime.values[1].degree;
    init_pol(&x_pow_l , l);
    set_coeff_pol(x_pow_l , l , 1);

    init_pol(&b , m);
    init_pol(&c , m);

    euclide_div_pol(&b , &r , BCprime.values[0] , x_pow_l);
    euclide_div_pol(&c , &r , BCprime.values[1] , x_pow_l);

    halfGCD(&Msecond , b , c);

    mul_matrix_pol(&temp , Mprime , M);
    mul_matrix_pol(Mgcd , Msecond , temp);

    destroy_pol(x_pow_m);
    destroy_pol(f);
    destroy_pol(g);
    destroy_pol(r);
    destroy_pol(Q);
    destroy_pol(Qneg);
    destroy_pol(x_pow_l);
    destroy_pol(b);
    destroy_pol(c);

    destroy_matrix_pol(M);
    destroy_matrix_pol(AB);
    destroy_matrix_pol(ABprime);
    destroy_matrix_pol(Mprime);
    destroy_matrix_pol(BCprime);
    destroy_matrix_pol(M_ABprime);
    destroy_matrix_pol(Msecond);
    destroy_matrix_pol(temp);
}
