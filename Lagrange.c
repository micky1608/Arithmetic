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
}

/* ************************************************************************************************************* */

void lagrange_interpolation(pol_bigQ *res , pol_bigQ_value *values , unsigned int nb_values) {

    change_degre_pol_bigQ(res , nb_values-1);
    set_all_coeffs_to_pol_bigQ_si(*res , 0 ,1);

    for(int j=0 ; j<nb_values ; j++) {

    }
}