//
// Created by root on 25/12/18.
//

#include "FFT.h"

void split_pol_bigfloat(pol_bigfloat *pol_even , pol_bigfloat *pol_odd , pol_bigfloat pol) {
    double degree = pol.degree;

    unsigned int degree_even , degree_odd;
    degree_odd = (unsigned int)floor(degree / 2);
    degree_even = pol.degree - degree_odd;

    change_degre_pol_bigfloat(pol_odd , degree_odd);
    change_degre_pol_bigfloat(pol_even , degree_even);

    unsigned int index_even = 0 , index_odd = 0;

    for(unsigned int i = 0 ; i<=pol.degree ; i++) {
        if(i%2 == 0) {
            mpfr_set(pol_even->coeffs[index_even] , pol.coeffs[i] , MPC_RNDDN);
            index_even++;
        }
        else {
            mpfr_set(pol_odd->coeffs[index_odd] , pol.coeffs[i] , MPC_RNDDN);
            index_odd++;
        }
    }
}

/* ********************************************************************************************************************** */

/**
 * Init and Fill the array the n-ieme roots of the unity
 * @param X
 * @param n
 */
void n_unit_roots(mpc_t X[] , unsigned int n) {
    for(unsigned int k=0 ; k<n ; k++) {
        mpc_init2(X[k] , MPC_RNDDN);
        mpc_set_d_d(X[k] , cos((2*k*M_PI)/n) , sin((2*k*M_PI)/n) , MPC_RNDDN);
    }


}