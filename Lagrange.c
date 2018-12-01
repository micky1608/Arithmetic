//
// Created by root on 30/11/18.
//

#include "Lagrange.h"

void init_polbigQ_value(pol_bigQ_value *polBigQValue , mpq_t Xi , mpq_t Yi) {
    mpq_init(polBigQValue->Xi);
    mpq_init(polBigQValue->Yi);

    mpq_set(polBigQValue->Xi , Xi);
    mpq_set(polBigQValue->Yi , Yi);
}

/* ************************************************************************************************************* */

void init_polbigQ_value_si(pol_bigQ_value *polBigQValue , int Xi_num , unsigned int Xi_den, int Yi_num , unsigned int Yi_den) {
    mpq_init(polBigQValue->Xi);
    mpq_init(polBigQValue->Yi);

    mpq_set_si(polBigQValue->Xi , Xi_num , Xi_den);
    mpq_set_si(polBigQValue->Yi , Yi_num , Yi_den);
}

/* ************************************************************************************************************* */

void destroy_polbigQ_value(pol_bigQ_value polBigQValue) {
    mpq_clears(polBigQValue.Xi , polBigQValue.Yi , NULL);
}

/* ************************************************************************************************************* */

void print_polbigQ_value(pol_bigQ_value polValue) {
    gmp_printf("Xi = %Qd\tYi = %Qd\n",polValue.Xi,polValue.Yi);
}

/* ************************************************************************************************************* */

void print_lagrange_constaints(pol_bigQ_value *constraints , unsigned int nb_constraint) {
    printf("\tLagrange constraints :\n");
    for(int i=0;i<nb_constraint;i++)
        print_polbigQ_value(constraints[i]);
    printf("\n");
}

/* ************************************************************************************************************* */

void lagrange_interpolation(pol_bigQ *res , pol_bigQ_value *values , unsigned int nb_values) {

    change_degre_pol_bigQ(res , nb_values-1);
    set_all_coeffs_to_pol_bigQ_si(*res , 0 ,1);


    pol_bigQ P , X_Xi , temp , resI , resI_num , resI_num2 , R;

    mpq_t Ui , Uij;
    mpq_inits(Ui , Uij , NULL);

    init_pol_bigQ(&P , 1);
    init_pol_bigQ(&X_Xi , 1);
    init_pol_bigQ(&temp , nb_values-1);
    init_pol_bigQ(&resI , nb_values-1);
    init_pol_bigQ(&resI_num , nb_values-1);
    init_pol_bigQ(&resI_num2 , nb_values-1);
    init_pol_bigQ(&R , nb_values-1);

    mpq_neg(P.coeffs[0] , values[0].Xi);
    mpq_set_d(P.coeffs[1] , 1);

    for(size_t i=1 ; i<nb_values ; i++) {
        mpq_neg(X_Xi.coeffs[0] , values[i].Xi);
        mpq_set_d(X_Xi.coeffs[1] , 1);

        mult_pol_bigQ(&temp , P , X_Xi);
        copy_pol_bigQ(&P , temp);
    }

    set_all_coeffs_to_pol_bigQ_si(temp , 0 , 1);

    for(size_t i=0 ; i<nb_values ; i++) {
        mpq_mult_polbigQ(&resI_num , P , values[i].Yi);

        mpq_neg(X_Xi.coeffs[0] , values[i].Xi);
        mpq_set_d(X_Xi.coeffs[1] , 1);

        mpq_set_d(Ui , 1);
        for(size_t j=0 ; j<nb_values ; j++) {
            if(j != i) {
                mpq_sub(Uij , values[i].Xi , values[j].Xi);
                mpq_mul(Ui , Ui , Uij);
            }
        }
        mpq_inv(Ui , Ui);

        mpq_mult_polbigQ(&resI_num2 , resI_num , Ui);

        euclideDiv_pol_bigQ(&resI , &R , resI_num2 , X_Xi);

        if(is_null_polbigQ(R) == FALSE)  {
            perror("Error in lagrange : Remainder not null !!!");
            exit(-1);
        }

        add_pol_bigQ(&temp , *res , resI);
        copy_pol_bigQ(res , temp);
    }


    mpq_clears(Ui , Uij , NULL);
    destroy_pol_bigQ(P);
    destroy_pol_bigQ(X_Xi);
    destroy_pol_bigQ(temp);
    destroy_pol_bigQ(resI);
    destroy_pol_bigQ(resI_num);
    destroy_pol_bigQ(resI_num2);
    destroy_pol_bigQ(R);
}