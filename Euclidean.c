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

void halfGCD(matrix_pol_double *Mgcd , pol_double A , pol_double B) {
    static int count = 0;
    printf("\t\t\t***** Half gcd call : %d *****\n",++count);

    print_pol_double(A , "Half gcd A");
    print_pol_double(B , "Half gcd B");
/*
    if(B.degree >= A.degree) {
        perror("HalfGCD B degree must be smaller");
        exit(EXIT_FAILURE);
    }
*/
    int n = A.degree , m = (int)ceil((double)n/2);

    printf("n : %d\tm : %d\n",n,m);
    
    if(n==0 || B.degree < m) {
        identity_matrix_pol_double(Mgcd , 2);
        return;
    }

    init_matrix_pol_double(Mgcd , 2 , 2);

    pol_double x_pow_m , f , g , r , Q , Qneg , x_pow_l , b , c;
    matrix_pol_double M , AB , ABprime , Mprime , BCprime , M_ABprime , Msecond , temp;

    init_pol_double(&x_pow_m , (unsigned int)m);
    set_coeff_pol_double(x_pow_m , (unsigned int)m , 1);

    init_pol_double(&f , m);
    init_pol_double(&g , m);
    init_pol_double(&r , m);
    init_pol_double(&Q , A.degree - B.degree);
    init_pol_double(&Qneg , A.degree-B.degree);

    init_matrix_pol_double(&AB , 2 , 1);
    init_matrix_pol_double(&ABprime , 2 , 1);
    init_matrix_pol_double(&Mprime , 2 , 1);
    init_matrix_pol_double(&BCprime , 2 , 1);
    init_matrix_pol_double(&M_ABprime , 2 , 2);
    init_matrix_pol_double(&temp , 2 , 2);

    setCoeff_matrix_pol_double(&AB , 0 , 0 , A);
    setCoeff_matrix_pol_double(&AB , 1 , 0 , B);

    euclide_div_pol_double(&f , &r , A , x_pow_m);
    euclide_div_pol_double(&g , &r , B , x_pow_m);

    print_pol_double(f , "f");
    print_pol_double(g , "g");

    halfGCD(&M , f , g); // recursive call

 //   print_matrix_pol_double(M , "M");   

    mul_matrix_pol_double(&ABprime , M , AB);

    if(ABprime.values[1].degree < m) {
        copy_matrix_pol_double(Mgcd , M);
        return;
    }
/*
    print_pol_double(ABprime.values[0] , "A'");
    print_pol_double(ABprime.values[1] , "B'");
*/
    euclide_div_pol_double(&Q , &r , ABprime.values[0] , ABprime.values[1]);

 //   print_pol_double(Q , "Q");

    scalar_mult_pol_double(&Qneg , Q , -1);
    set_coeff_constant_matrix_pol_double(&M_ABprime , 0 , 1 , 1);
    set_coeff_constant_matrix_pol_double(&M_ABprime , 1 , 0 , 1);
    setCoeff_matrix_pol_double(&M_ABprime , 1 , 1 , Qneg);

    mul_matrix_pol_double(&BCprime , M_ABprime , ABprime);

/*
    print_matrix_pol_double(M_ABprime , "M_AB'");
    print_matrix_pol_double(BCprime , "BC'");
*/
    /* ----------------------------------------------------------------- */

    int l = 2*m - ABprime.values[1].degree;
    init_pol_double(&x_pow_l , l);
    set_coeff_pol_double(x_pow_l , l , 1);

    init_pol_double(&b , m);
    init_pol_double(&c , m);
/*
    print_pol_double(BCprime.values[0] , "B'");
    print_pol_double(BCprime.values[1] , "C'");
    print_pol_double(x_pow_l , "X^l");
*/
    euclide_div_pol_double(&b , &r , BCprime.values[0] , x_pow_l);
    euclide_div_pol_double(&c , &r , BCprime.values[1] , x_pow_l);
/*
    print_pol_double(b , "b");
    print_pol_double(c , "c");
*/
    halfGCD(&Msecond , b , c); // recursive call

/*
    print_matrix_pol_double(Msecond , "M''");
    print_matrix_pol_double(Mprime , "M'");
    print_matrix_pol_double(M , "M");
*/

    mul_matrix_pol_double(&temp , M_ABprime , M);
    mul_matrix_pol_double(Mgcd , Msecond , temp);

    destroy_pol_double(x_pow_m);
    destroy_pol_double(f);
    destroy_pol_double(g);
    destroy_pol_double(r);
    destroy_pol_double(Q);
    destroy_pol_double(Qneg);
    destroy_pol_double(x_pow_l);
    destroy_pol_double(b);
    destroy_pol_double(c);

    destroy_matrix_pol_double(M);
    destroy_matrix_pol_double(AB);
    destroy_matrix_pol_double(ABprime);
    destroy_matrix_pol_double(Mprime);
    destroy_matrix_pol_double(BCprime);
    destroy_matrix_pol_double(M_ABprime);
    destroy_matrix_pol_double(Msecond);
    destroy_matrix_pol_double(temp);
}

/* ********************************************************************************************************************** */

void fast_euclide(matrix_pol_double *Mab , pol_double A , pol_double B) {

    static int count = 0;
    printf("\t***** Fast euclide call : %d *****\n",++count);

    pol_double Q , r;
    matrix_pol_double Mgcd , AB , R , M , R2 , MR, temp;
    
    init_matrix_pol_double(Mab , 2 , 2);
    init_matrix_pol_double(&AB , 2 , 1);
    init_matrix_pol_double(&R , 2 , 1);
   
    halfGCD(&Mgcd , A , B);

    print_matrix_pol_double(Mgcd , "Mgcd");

    setCoeff_matrix_pol_double(&AB , 0 , 0 , A);
    setCoeff_matrix_pol_double(&AB , 1 , 0 , B);

    mul_matrix_pol_double(&R , Mgcd , AB);

    if(is_zero_pol_double(R.values[1])) { 
        copy_matrix_pol_double(Mab , Mgcd);
        return;
    }

    init_matrix_pol_double(&M , 2 , 2);
    init_matrix_pol_double(&R2 , 2 , 1);

    init_pol_double(&Q , 0);
    init_pol_double(&r , 0);

    set_coeff_constant_matrix_pol_double(&M , 0 , 0 , 0);
    set_coeff_constant_matrix_pol_double(&M , 0 , 1 , 1);
    set_coeff_constant_matrix_pol_double(&M , 1 , 0 , 1);

    print_matrix_pol_double(R , "R");
/*
    print_pol_double(R.values[0] , "R[0,0]");
    print_pol_double(R.values[1] , "R[1,0]");
*/  
    euclide_div_pol_double(&Q , &r , R.values[0] , R.values[1]);
    print_pol_double(Q , "Q");
    scalar_mult_pol_double(&M.values[3] , Q , -1);
    
    print_matrix_pol_double(M , "M");
    mul_matrix_pol_double(&R2 , M , R);

    //if(is_zero_pol_double(R2.values[1])) {
    if(R2.values[1].degree == 0) {
        mul_matrix_pol_double(Mab , M , Mgcd);
        return;
    }

    print_matrix_pol_double(R2 , "R2");

    fast_euclide(&MR , R2.values[0] , R2.values[1]); // recursive call

    init_matrix_pol_double(&temp , 2 , 2);

    mul_matrix_pol_double(&temp , M , Mgcd);
    mul_matrix_pol_double(Mab , MR , temp);

    
    destroy_pol_double(Q);
    destroy_pol_double(r);

    destroy_matrix_pol_double(Mgcd);
    destroy_matrix_pol_double(AB);
    destroy_matrix_pol_double(R);
    destroy_matrix_pol_double(M);
    destroy_matrix_pol_double(R2);
    destroy_matrix_pol_double(MR);
    destroy_matrix_pol_double(temp);
}
