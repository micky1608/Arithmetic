//
// Created by root on 25/12/18.
//

#include "FFT.h"
#include "Complex.h"

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

/* ********************************************************************************************************************** */

void polynomial_naive_evaluation(mpc_t *y , pol_bigfloat P , mpc_t x) {
    mpc_set_d(*y , 0 , MPC_RNDDN);

    mpc_t temp;
    mpc_init2(temp , 64);
    for(unsigned int k=0 ; k<= P.degree ; k++) {
        mpc_pow_d(temp , x , k , MPC_RNDDN);
        mpc_mul_fr(temp , temp , P.coeffs[k] , MPC_RNDDN);
        mpc_add(*y , *y , temp , MPC_RNDDN);
    }
    mpc_clear(temp);
}

/* ********************************************************************************************************************** */

/**
 * Evaluate the polynomial P in the set on points Xi.
 * This set contains the n-ieme roots of unity. The second half of the set contains the negatives roots of the first half
 * @param P_eval
 * @param P
 * @param Xi
 * @param n
 */
void FFT_evaluation(mpc_t *P_eval , pol_bigfloat P , mpc_t *Xi , unsigned int n) {
    if(log_base_2(n) - floor(log_base_2(n)) != 0) {
        perror("FFT_evaluation : n must be a power of 2");
        return;
    }

    if(n == 1) {
        polynomial_naive_evaluation(&P_eval[0] , P , Xi[0]);
        return;
    }

    mpc_t Xi2[n/2];
    mpc_t eval_even[n/2] , eval_odd[n/2];

    for(unsigned int k=0 ; k<n/2 ; k++) {
        mpc_init2(Xi2[k] , 64);
        mpc_pow_d(Xi2[k] , Xi[k] , 2 , MPC_RNDDN);

        mpc_init2(eval_even[k] , MPC_RNDDN);
        mpc_init2(eval_odd[k] , MPC_RNDDN);
    }

    pol_bigfloat P_even,P_odd;
    init_pol_bigfloat(&P_even , P.degree / 2);
    init_pol_bigfloat(&P_odd , P.degree / 2);

    split_pol_bigfloat(&P_even , &P_odd , P);
    
    FFT_evaluation(eval_even , P_even , Xi2 , n/2);
    FFT_evaluation(eval_odd , P_odd , Xi2 , n/2);

    for(unsigned int k=0 ; k<n/2 ; k++) {
        // eval[k] = eval_even[k] + Xi[k]*eval_odd[k]
        // eval[n/2 + k] = eval_even[k] - Xi[k]*eval_odd[k]

        mpc_init2(P_eval[k] , MPC_RNDDN);
        mpc_init2(P_eval[n/2 + k] , MPC_RNDDN);

        mpc_mul(P_eval[k] , eval_odd[k] , Xi[k] , MPC_RNDDN);
        mpc_add(P_eval[k] , P_eval[k] , eval_even[k] , MPC_RNDDN);

        mpc_mul(P_eval[n/2 + k] , eval_odd[k] , Xi[k] , MPC_RNDDN);
        mpc_neg(P_eval[n/2 + k] , P_eval[n/2 + k] , MPC_RNDDN);
        mpc_add(P_eval[n/2 + k] , P_eval[n/2 + k] , eval_even[k] , MPC_RNDDN);
    }

    destroy_pol_bigfloat(P_even);
    destroy_pol_bigfloat(P_odd);
}

/* ********************************************************************************************************************** */

void FFT_multiplication(mpc_t *R_eval , mpc_t *P_eval , mpc_t Q_eval) {
    //TODO
}

/* ********************************************************************************************************************** */

void FFT_reverse(pol_bigfloat *R , mpc_t *R_eval , mpc_t omega , unsigned int n) {
    //TODO
}

/* ********************************************************************************************************************** */

void mul_pol_bigfloat_FFT(pol_bigfloat *R , pol_bigfloat P , pol_bigfloat Q) {
    //TODO
}