/*
* Created on 2019-05-16
*/

#ifndef HEADER_SERIE_H
#define HEADER_SERIE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "arithmetic.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


/*
    1 : long
    2 : double
*/
#define S_COEFF_T 1

/* *************************************************************************************************************************************************** */

#if S_COEFF_T == 1
    #define S_COEFF_LONG
    typedef long s_coeff_t;

#elif S_COEFF_T == 2
    #define S_COEFF_DOUBLE
    typedef double s_coeff_t;
#endif

/* *************************************************************************************************************************************************** */


typedef struct serie {
    unsigned int precision;
    s_coeff_t *coeffs;
} serie;    

void init_serie(serie *s, unsigned int precision);

void set_all_coeffs_serie(serie *s, s_coeff_t *coeffs, size_t size_coeffs_array);

void destroy_serie(serie s);

void change_precision_serie(serie *s, unsigned int new_precision);

void print_serie(serie s , char *name);

void reduce_serie_ff(serie s , long p);

void add_serie(serie *s , serie a , serie b);
void add_serie_ff(serie *s , serie a , serie b , long p);

void mul_serie(serie *s , serie a , serie b);
void mul_serie_ff(serie *s , serie a , serie b , long p);

void inv_serie_ff(serie *inv , serie *s , unsigned int precision , long p);

#endif