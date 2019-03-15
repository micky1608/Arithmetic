#include "pol_double.h"

/* ********************************************************************************************************************** */

/**
 * Allocate the memory for the coefficients
 * @param degre
 */
void init_pol_double(pol_double *polynomial, unsigned int degree_p) {
    if(degree_p >=0) {
        polynomial->degree = degree_p;
        polynomial->coeffs = malloc((degree_p + 1) * sizeof(double));
        memset(polynomial->coeffs , 0 , (degree_p+1)*sizeof(double));
    }
    else {
        perror("Try to create a polynomial with negative degree !!");
        exit(-1);
    }
}

/* ********************************************************************************************************************** */

/**
 * Free the memory allocated for a polynomial
 * @param polynomial
 */
void destroy_pol_double(pol_double polynomial) {
    free(polynomial.coeffs);
}

/* ********************************************************************************************************************** */

/**
 * Change the coefficients of a polynomial from an array containing [ f0, f1 .... ]
 * @param polynomial
 * @param coeffs
 */
void set_all_coeffs_pol_double(pol_double polynomial, double *coeffs, size_t size_coeffs_array) {
    for(int i=0 ; i<size_coeffs_array ; i++) {
        if (i <= polynomial.degree)
            polynomial.coeffs[i] = coeffs[i];
        else
            break;
    }
}

/* ********************************************************************************************************************** */

/**
 * Define all the coefficients of a polynomial between -max and max
 * @param polynomial
 * @param max
 */
void set_all_coeffs_random_pol_double(pol_double polynomial, unsigned int max) {

    for(int i=0 ; i<=polynomial.degree ; i++) {
        int fi_i = (rand() % (2*max)) - max;
        while(i==polynomial.degree && fi_i==0) {
            fi_i = (rand() % (2*max)) - max;
        }
        polynomial.coeffs[i] = fi_i;
    }
}

/* ********************************************************************************************************************** */

/**
 * Set all the coefficients to the same value
 * @param polynomial
 * @param value
 */
void set_all_coeffs_to_pol_double(pol_double polynomial, double value) {
    for(int i=0 ; i<=polynomial.degree ; i++)
        polynomial.coeffs[i] = value;

}

/* ********************************************************************************************************************** */

/**
 * Change one coefficient in a polynomial
 * @param polynome
 * @param degre_coeff
 * @param newValue
 */
void set_coeff_pol_double(pol_double polynomial, unsigned int degree_coeff, double newValue) {
    if(degree_coeff <= polynomial.degree)
        polynomial.coeffs[degree_coeff] = newValue;
}

/* ********************************************************************************************************************** */


/**
 * Change the degree of a polynomial
 * If the new degree is greater than the old, the function complete with 0
 * If the new degree is smaller than the old, some coefficients are lost
 * @param polynomial
 */
void change_degre_pol_double(pol_double *polynomial, unsigned int new_degree) {
    if(new_degree == polynomial->degree || new_degree < 0) return;

    double *new_coeffs = calloc(new_degree+1 , sizeof(double));

    if(new_degree < polynomial->degree) {
        for(int i=0 ; i<=polynomial->degree ; i++) {
            if(i<=new_degree)
                new_coeffs[i] = polynomial->coeffs[i];
        }
    }
    else {
            for(int i=0 ; i<=new_degree ; i++) {
                if(i<=polynomial->degree) {
                    new_coeffs[i] = polynomial->coeffs[i];
                }
                else
                    new_coeffs[i] = 0;
            }
    }

    free(polynomial->coeffs);
    polynomial->coeffs = new_coeffs;
    polynomial->degree = new_degree;
}

/* ********************************************************************************************************************** */

/**
 * Copy the polynomial src in res
 * @param res
 * @param polynomial
 */
void copy_pol_double(pol_double *res, pol_double polynomial) {
    change_degre_pol_double(res, polynomial.degree);
    for(int i=0 ; i<=res->degree ; i++)
        res->coeffs[i] = polynomial.coeffs[i];
}

/* ********************************************************************************************************************** */

void print_pol_double(pol_double polynomial, char *name) {
    if(name != NULL) printf("%s (degree %d): ",name,polynomial.degree);
    if(is_zero_pol_double(polynomial)) {
        printf("0");
        if(name != NULL) printf("\n");
        return;
    }

    double temp;

    if (polynomial.coeffs[0]) printf("%.3f",polynomial.coeffs[0]);
    for(int i=1 ; i <= polynomial.degree ; i++) {

        if (polynomial.coeffs[i] < 0 )
            printf(" - ");
        else if (polynomial.coeffs[i] > 0 )
            printf(" + ");
        else
            continue;

        temp = (polynomial.coeffs[i] > 0) ? polynomial.coeffs[i] : -1*polynomial.coeffs[i];
        printf("(");
        printf("%.3f",temp);
        printf(" * X^%d)" , i);

    }
    
    if(name != NULL) printf("\n");
}

/* ********************************************************************************************************************** */

void print_pol_double_center(pol_double polynomial , char *name , unsigned int size) {
    if(polynomial.degree >= size) {
        print_pol_double(polynomial , name);
        return;
    }

    int blank , i=0;

    blank = size - polynomial.degree; 

    while(i<ceil((double)blank/2)) {
        if(!i) printf("   ");
        else printf("         ");
        i++;
    }
    print_pol_double(polynomial , name);
    while(i<=blank) {
        printf("         ");
        i++;
    }

}

/* ********************************************************************************************************************** */

/**
 * Add two polynomials
 * @param res
 * @param A
 * @param B
 */
void add_pol_double(pol_double *res, pol_double A, pol_double B) {

    if(res->degree != MAX(A.degree , B.degree))
        change_degre_pol_double(res, MAX(A.degree, B.degree));

    for(int i=0 ; i<=res->degree ; i++) {
        if(i <= A.degree) {
            if(i <= B.degree) {
               res->coeffs[i] = A.coeffs[i] + B.coeffs[i];
            }
            else
                set_coeff_pol_double(*res, i, A.coeffs[i]);
        }
        else
            set_coeff_pol_double(*res, i, B.coeffs[i]);
    }
}

/* ********************************************************************************************************************** */

/**
 * Subtract two polynomials
 * @param res
 * @param A
 * @param B
 */
void sub_pol_double(pol_double *res, pol_double A, pol_double B) {

    if(res->degree != MAX(A.degree , B.degree))
        change_degre_pol_double(res, MAX(A.degree, B.degree));

    for(int i=0 ; i<=res->degree ; i++) {
        if(i <= A.degree) {
            if(i <= B.degree) {
                res->coeffs[i] = A.coeffs[i] - B.coeffs[i];
            }
            else
                set_coeff_pol_double(*res, i, A.coeffs[i]);
        }
        else
            *(res->coeffs+i) = -1 * B.coeffs[i];
    }

}

/* ********************************************************************************************************************** */

/**
 * Naive multiplication of two polynomials
 * @param res
 * @param A
 * @param B
 */
void mult_pol_double(pol_double *res, pol_double A, pol_double B) {
    if(is_zero_pol_double(A) || is_zero_pol_double(B)) {
        change_degre_pol_double(res , 0);
        set_coeff_pol_double(*res , 0 , 0);
        return;
    }
    change_degre_pol_double(res, A.degree + B.degree);
    double temp;

    for(int i=0 ; i<=res->degree ; i++) {
        res->coeffs[i] = 0;

        int index_A = 0;
        int index_B = i;

        while (index_A <= i && index_B >=0) {

            if(index_A <= A.degree && index_B <= B.degree) {
                temp = A.coeffs[index_A] * B.coeffs[index_B];
                res->coeffs[i] += temp;
            }
            index_A++;
            index_B--;
        }
    }
}

/* ********************************************************************************************************************** */

void scalar_mult_pol_double(pol_double *res , pol_double A , double lambda) {
    change_degre_pol_double(res , A.degree);
    for(unsigned int i=0 ; i<=res->degree ; i++)
        set_coeff_pol_double(*res , i , lambda*A.coeffs[i]);
}

/* ********************************************************************************************************************** */


int is_zero_pol_double(pol_double polynomial) {
    for(unsigned int i=0 ; i<=polynomial.degree ; i++)
        if (polynomial.coeffs[i]) return 0;
    return 1;
}

/* ********************************************************************************************************************** */


void euclide_div_pol_double(pol_double *Q , pol_double *R , pol_double A , pol_double B) {

   if(A.degree < B.degree) {
        perror("Euclidean division : A.degree < B.degree !!");
        return;
    }

    if(is_zero_pol_double(B)) {
        perror("Euclidean division : B can't be 0 !!");
        return;
    }

    double *coeffsA = (double*) malloc(sizeof(double) * (A.degree+1));
    double *coeffsB = (double*) malloc(sizeof(double) * (B.degree+1));
    double *coeffsQ = (double*) malloc(sizeof(double) * (A.degree - B.degree +1));
    double *coeffsR = (double*) malloc(sizeof(double) * (A.degree+1));
    double *coeffsT = (double*) malloc(sizeof(double) * (A.degree+1));

    int sizeB = B.degree, 
        sizeR = A.degree, 
        sizeQ = A.degree - B.degree;

    for(int i=0 ; i<=A.degree ; i++) {
        coeffsA[i] = (double)A.coeffs[i];
        coeffsR[i] = (double)A.coeffs[i];
        coeffsT[i] = 0;
        if(i <= sizeB) coeffsB[i] = (double)B.coeffs[i];
        if(i <= sizeQ) coeffsQ[i] = 0;
    }

    double a,b,t;

    b = coeffsB[sizeB];

    while(sizeR >= sizeB) {
        a = coeffsR[sizeR];
        t = a/b;

        coeffsQ[sizeR - sizeB] += t;

        for(int i=0 ; i<=sizeR ; i++) {
            if(i < (sizeR-sizeB)) coeffsT[i] = 0;
            else coeffsT[i] = t*coeffsB[i - (sizeR - sizeB)];
        }

        for(int i=0 ; i<=sizeR ; i++) coeffsR[i] -= coeffsT[i];
        for(;coeffsR[sizeR] == 0 && sizeR>0; --sizeR);

    }

    for(;coeffsQ[sizeQ] == 0 && sizeQ>0; --sizeQ);

    change_degre_pol_double(Q , sizeQ);
    change_degre_pol_double(R , sizeR);
    
    for(int i=0 ; i<=MAX(sizeR,sizeQ) ; i++) {
        if(i<=sizeR) R->coeffs[i] = (double)coeffsR[i];
        if(i<=sizeQ) Q->coeffs[i] = (double)coeffsQ[i];
    }

    free(coeffsA);
    free(coeffsB);
    free(coeffsT);
    free(coeffsR);
    free(coeffsQ);
}

/* ********************************************************************************************************************** */


void horner_eval_double(double *res , pol_double f , double x) {
    *res = f.coeffs[f.degree];
    for(unsigned int i=1 ; i<=f.degree ; i++) {
        *res *= x;
        *res += f.coeffs[f.degree - i];
    }
}

void horner_eval_multi_double(double *res , pol_double f , double *x , unsigned int nb_x) {
    for(unsigned int i=0 ; i<nb_x ; i++) horner_eval_double(res + i , f , *(x + i));
}

