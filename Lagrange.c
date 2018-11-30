//
// Created by root on 30/11/18.
//

#include "Lagrange.h"

void lagrange_interpolation(pol_bigQ *res , pol_value *values , unsigned int nb_values) {

    change_degre_pol_bigQ(res , nb_values-1);
    set_all_coeffs_to_pol_bigQ_si(*res , 0 ,1);

    for(int j=0 ; j<nb_values ; j++) {

    }
}