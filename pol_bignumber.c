//
// Created by root on 29/10/18.
//

#include "pol_bignumber.h"

/**
 * Copy a polynomial to create a float equivalent
 * @param res
 * @param A
 */
void pol_bigint_to_pol_bigfloat(pol_bigfloat *res , pol_bigint A) {
    if(res->degree != A.degree)
        change_degre_pol_bigfloat(res , A.degree);

    for(int i=0 ; i<=res->degree ; i++) {
        mpf_set_z(res->coeffs[i] , A.coeffs[i]);
        gmp_printf("A.coeff[%d] = %Zd\tres->coeffs[%d] = %.2F\n",i,A.coeffs[i],i,res->coeffs[i]);
    }
}

/* ********************************************************************************************************************** */

/**
 * Copy a polynomial to create an int equivalent.
 * The coefficients are truncated.
 * @param res
 * @param A
 */
void pol_bigfloat_to_pol_bigint(pol_bigint *res , pol_bigfloat A) {
    if(res->degree != A.degree)
        change_degre_pol_bigint(res , A.degree);

    for(int i=0 ; i<=res->degree ; i++)
        mpz_set_f(res->coeffs[i] , A.coeffs[i]);
}
