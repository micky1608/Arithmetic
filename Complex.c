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