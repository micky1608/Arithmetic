//
// Created by root on 24/12/18.
//


#include <mpfr.h>
#include "Complex.h"

/**
 * Print a complex number with the given precision.
 * If the precision is NULL, then the number is displayed with the NN precision.
 * @param precision
 */
void print_complex(mpc_t z , char *name) {
    mpfr_t realPart , ImPart;
    mpfr_inits(realPart,ImPart,NULL);

    mpc_real(realPart , z , MPC_RNDDN);
    mpc_imag(ImPart , z , MPC_RNDDN);

    printf("%s : ",name);

    mpfr_printf("%.3Rf " , realPart);

    if(mpfr_cmp_d(ImPart , 0) < 0) printf("- i * ");
    else printf("+ i * ");

    mpfr_abs(ImPart , ImPart , MPC_RNDDN);

    mpfr_printf("%.3Rf\n" , ImPart);

    mpfr_clears(realPart,ImPart,NULL);
}

/* ********************************************************************************************************************** */

void init_pol_complex(pol_complex *pol , unsigned int degree) {
    pol->degree = degree;
    pol->coeffs = calloc(degree + 1 , sizeof(mpc_t));
    for(int i=0 ; i<= pol->degree ; i++) mpc_init2(pol->coeffs[i] , 64);
}

/* ********************************************************************************************************************** */

void set_coeffs_pol_complex(pol_complex *pol , mpc_t *new_coeffs , unsigned int nb_coeffs) {
    if(nb_coeffs > pol->degree) return;
    for(unsigned int i=0 ; i<=nb_coeffs ; i++) {
        mpc_set(pol->coeffs[i] , new_coeffs[i] , MPC_RNDDN);
    }
}

/* ********************************************************************************************************************** */

void change_degre_pol_complex(pol_complex *pol, unsigned int new_degree) {
    if(new_degree == pol->degree || new_degree < 0) return;

    mpc_t *new_coeffs = calloc(new_degree+1 , sizeof(mpc_t));

    if(new_degree < pol->degree) {
        for(int i=0 ; i<=pol->degree ; i++) {
            if(i<=new_degree) {
                mpc_init2(new_coeffs[i] , 64);
                mpc_set(new_coeffs[i] , pol->coeffs[i] , MPC_RNDDN);
            }
            mpc_clear(pol->coeffs[i]);
        }
    }
    else {
        for(int i=0 ; i<=new_degree ; i++) {
            if(i<=pol->degree) {
                mpc_init2(new_coeffs[i] , 64);
                mpc_set(new_coeffs[i] , pol->coeffs[i] , MPC_RNDDN);
                mpc_clear(pol->coeffs[i]);
            }
            else {
                mpc_init2(new_coeffs[i] , 64);
                mpc_set_d(new_coeffs[i] , 0 , MPC_RNDDN);
            }
        }
    }

    free(pol->coeffs);
    pol->coeffs = new_coeffs;
    pol->degree = new_degree;
}

/* ********************************************************************************************************************** */

void destroy_pol_complex(pol_complex pol) {
    for(int i=0 ; i<= pol.degree ; i++)
        mpc_clear(pol.coeffs[i]);
    free(pol.coeffs);
}