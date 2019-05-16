//
// Created by micky on 27/02/19.
//

#include "test.h"

/* ********************************************************************************************************************** */


void test_pol() {
    pol_double A,B,Q,R,F;

    //init_pol_double(&A , 3);
    init_pol_double(&A , 2);

    init_pol_double(&B , 1);
    init_pol_double(&Q , 1);
    init_pol_double(&R , 1);
    init_pol_double(&F , 0);

    //double coeffA[] = {-4,0,-2,1};
    //double coeffB[] = {-3,1};

    double coeffA[] = {2, -1 , 3};
    double coeffB[] = {6 , 4};

    //set_all_coeffs_pol_double(A , coeffA , 4);
    set_all_coeffs_pol_double(A , coeffA , 3);
    set_all_coeffs_pol_double(B , coeffB , 2);

    euclide_div_pol_double(&Q , &R , A , B);

    print_pol_double(A , "A");
    print_pol_double(B , "B");
    print_pol_double(Q , "Q");
    print_pol_double(R , "R");

    copy_pol_double(&F , A);
    print_pol_double(F , "F (copy A)");
    
    mult_pol_double(&F , A , B);
    print_pol_double(F , "F (A*B)");

    add_pol_double(&F , Q , B);
    print_pol_double(F , "F (Q+B)");


    destroy_pol_double(A);
    destroy_pol_double(B);
    destroy_pol_double(Q);
    destroy_pol_double(R);
    destroy_pol_double(F);

    pol D, D_derivate;
    init_pol(&D , 5);
    long coeffsD[] = {1,2,3,4,5,6};
    set_all_coeffs_pol(D , coeffsD , 6);

    derivate_pol(&D_derivate , D);

    print_pol(D , "D");
    print_pol(D_derivate , "D'");

    reduce_pol_ff(&D , 3);

    print_pol(D , "D in F3[x]");

    destroy_pol(D);
    destroy_pol(D_derivate);

    pol E,G,Q2,R2;
    init_pol(&E , 3);
    init_pol(&G , 2);

    long coeffsE[] = {1,1,1,2};
    long coeffsF[] = {1,0,2};

    set_all_coeffs_pol(E , coeffsE , 4);
    set_all_coeffs_pol(G , coeffsF , 3);


    euclide_div_pol_ff(&Q2 , &R2 , E , G , 3);

    print_pol(Q2 , "Q2");
    print_pol(R2 , "R2");


    destroy_pol(E);
    destroy_pol(G);
    destroy_pol(Q2);
    destroy_pol(R2);


    pol X,Y,gcd;
    init_pol(&X , 2);
    init_pol(&Y , 2);

    long coeffsX[] = {0,1,2};
    long coeffsY[] = {1,0,2};

    set_all_coeffs_pol(X , coeffsX , 3);
    set_all_coeffs_pol(Y , coeffsY , 3);

    gcd_pol_ff(&gcd , X , Y , 3);

    print_pol(X , "X in F3[x]");
    print_pol(Y , "Y in F3[x]");
    print_pol(gcd , "gcd(X,Y) in F3[x]");


    destroy_pol(X);
    destroy_pol(Y);
    destroy_pol(gcd);



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
    matrix_pol_double Mgcd;
    pol_double A,B;
    
    init_pol_double(&A , 2);
    init_pol_double(&B , 1);

    double coeffsA[] = {2, -1 , 3};
    double coeffsB[] = {6 , 4};

    set_all_coeffs_pol_double(A , coeffsA , 3);
    set_all_coeffs_pol_double(B , coeffsB , 2);
    
    halfGCD(&Mgcd , A , B);

    print_matrix_pol_double(Mgcd , "Mgcd");

    destroy_matrix_pol_double(Mgcd);
    destroy_pol_double(A);
    destroy_pol_double(B);
}

/* ********************************************************************************************************************** */

void test_fastEuclide() {
    pol_double A,B,Q,R;

    matrix_pol_double Mab, AB, Rm;

    //init_pol_double(&A , 3);
    init_pol_double(&A , 2);

    init_pol_double(&B , 1);
    init_pol_double(&Q , 1);
    init_pol_double(&R , 1);

    init_matrix_pol_double(&AB , 2 , 1);
    init_matrix_pol_double(&Rm , 2 , 1);

    //double coeffA[] = {-4,0,-2,1};
    //double coeffB[] = {-3,1};

    double coeffA[] = {2, -1 , 3};
    double coeffB[] = {6 , 4};

    //set_all_coeffs_pol_double(A , coeffA , 4);
    set_all_coeffs_pol_double(A , coeffA , 3);
    set_all_coeffs_pol_double(B , coeffB , 2);

    fast_euclide(&Mab , A , B);

    print_pol_double(A , "A");
    print_pol_double(B , "B");
    
    print_matrix_pol_double(Mab , "Mab");

    setCoeff_matrix_pol_double(&AB , 0 , 0 , A);
    setCoeff_matrix_pol_double(&AB , 1 , 0 , B);

    mul_matrix_pol_double(&Rm , Mab , AB);

    print_matrix_pol_double(Rm , "Rm");

    destroy_pol_double(A);
    destroy_pol_double(B);
    destroy_pol_double(Q);
    destroy_pol_double(R);

    destroy_matrix_pol_double(Mab);
    destroy_matrix_pol_double(Rm);
    destroy_matrix_pol_double(AB);
}

/* ********************************************************************************************************************** */



void test_pol_fact() {
    pol F;
    init_pol(&F , 10);

    long coeffsA[] = {2,0,1,0,1,2,0,2,1,2,1};

    set_all_coeffs_pol(F , coeffsA , 11);

    fact_list_info *L;

    fact_algo2(&L , F , 3);

    print_fact_list(L);

    clean_fact_list(L);
}


void test_serie() {
    serie s;
    init_serie(&s , 6);
    s_coeff_t coeffs[] = {0,0,3,4,5,6};

    set_all_coeffs_serie(&s , coeffs , 6);
    
    print_serie(s , "S");

    destroy_serie(s);


}
