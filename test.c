//
// Created by micky on 27/02/19.
//

#include "test.h"

/* ********************************************************************************************************************** */

void test_pol() {
    pol A,B,Q,R,F;
    init_pol(&A , 3);
    init_pol(&B , 1);
    init_pol(&Q , 1);
    init_pol(&R , 1);
    init_pol(&F , 0);

    long coeffA[] = {-4,0,-2,1};
    long coeffB[] = {-3,1};

    set_all_coeffs_pol(A , coeffA , 4);
    set_all_coeffs_pol(B , coeffB , 2);

    euclide_div_pol(&Q , &R , A , B);

    print_pol(A , "A");
    print_pol(B , "B");
    print_pol(Q , "Q");
    print_pol(R , "R");

    copy_pol(&F , A);
    print_pol(F , "F (copy A)");
    
    mult_pol(&F , A , B);
    print_pol(F , "F (A*B)");

    add_pol(&F , Q , B);
    print_pol(F , "F (Q+B)");


    destroy_pol(A);
    destroy_pol(B);
    destroy_pol(Q);
    destroy_pol(R);
    destroy_pol(F);

}

/* ********************************************************************************************************************** */

void test_horner() {
    pol f;
    init_pol(&f , 3);

    long coeffs[] = {1,2,3,4};
    set_all_coeffs_pol(f , coeffs , 4);
    print_pol(f , "f");

    long res[3];
    long x[] = {1,2,-1};
    horner_eval_multi(res , f , x , 3);
    for(unsigned int i=0 ; i<3 ; i++)printf("f(%ld) = %ld\n",x[i] , res[i]);
    destroy_pol(f);
}

/* ********************************************************************************************************************** */

void test_matrix_pol() {

    matrix_pol M;
    init_matrix_pol(&M , 2 , 2);

    pol A;
    init_pol(&A , 2);
    long coeffs[] = {1,2,3};
    set_all_coeffs_pol(A , coeffs , 3);
    set_allCoeff_matrix_pol(&M , A);

    print_matrix_pol(M , "M");

    change_degre_pol(&M.values[2] , 1);
    change_degre_pol(&M.values[1] , 3);
    set_coeff_pol(M.values[1] , 3 , 4);

    print_matrix_pol(M , "M");

    change_degre_pol(&M.values[0] , 0);
    set_coeff_pol(M.values[0] , 0 , 0);

    print_matrix_pol(M , "M");

    matrix_pol id;
    identity_matrix_pol(&id , 2);

    matrix_pol res;
    mul_matrix_pol(&res , M , M);

    print_matrix_pol(res , "M*M");

    destroy_pol(A);
    destroy_matrix_pol(M);
}

/* ********************************************************************************************************************** */

void test_halfGCD() {
    matrix_pol Mgcd;
    pol A,B;
    
    init_pol(&A , 2);
    init_pol(&B , 1);

    long coeffsA[] = {2, -1 , 3};
    long coeffsB[] = {6 , 4};

    set_all_coeffs_pol(A , coeffsA , 3);
    set_all_coeffs_pol(B , coeffsB , 2);
    
    halfGCD(&Mgcd , A , B);

    destroy_matrix_pol(Mgcd);
    destroy_pol(A);
    destroy_pol(B);
}
