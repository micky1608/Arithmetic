#include "pol.h"

/* ********************************************************************************************************************** */

/**
 * Allocate the memory for the coefficients
 * @param degre
 */
void init_pol(pol *polynomial, unsigned int degree_p) {
    if(degree_p >=0) {
        polynomial->degree = degree_p;
        polynomial->coeffs = malloc((degree_p + 1) * sizeof(long));
        memset(polynomial->coeffs , 0 , (degree_p+1)*sizeof(long));
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
void destroy_pol(pol polynomial) {
    free(polynomial.coeffs);
}

/* ********************************************************************************************************************** */

/**
 * Change the coefficients of a polynomial from an array containing [ f0, f1 .... ]
 * @param polynomial
 * @param coeffs
 */
void set_all_coeffs_pol(pol polynomial, long *coeffs, size_t size_coeffs_array) {
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
void set_all_coeffs_random_pol(pol polynomial, unsigned int max) {

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
void set_all_coeffs_to_pol(pol polynomial, long value) {
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
void set_coeff_pol(pol polynomial, unsigned int degree_coeff, long newValue) {
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
void change_degre_pol(pol *polynomial, unsigned int new_degree) {
    if(new_degree == polynomial->degree || new_degree < 0) return;

    long *new_coeffs = calloc(new_degree+1 , sizeof(long));

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
void copy_pol(pol *res, pol polynomial) {
    change_degre_pol(res, polynomial.degree);
    for(int i=0 ; i<=res->degree ; i++)
        res->coeffs[i] = polynomial.coeffs[i];
}

/* ********************************************************************************************************************** */

void print_pol(pol polynomial, char *name) {
    if(name != NULL) printf("%s (degree %d): ",name,polynomial.degree);
    if(is_zero_pol(polynomial)) {
        printf("0");
        if(name != NULL) printf("\n");
        return;
    }

    long temp;

    if (polynomial.coeffs[0]) printf("%ld",polynomial.coeffs[0]);
    for(int i=1 ; i <= polynomial.degree ; i++) {

        if (polynomial.coeffs[i] < 0 )
            printf(" - ");
        else if (polynomial.coeffs[i] > 0 )
            printf(" + ");
        else
            continue;

        temp = (polynomial.coeffs[i] > 0) ? polynomial.coeffs[i] : -1*polynomial.coeffs[i];
        printf("(");
        printf("%ld",temp);
        printf(" * X^%d)" , i);

    }
    
    if(name != NULL) printf("\n");
}

/* ********************************************************************************************************************** */

void print_pol_center(pol polynomial , char *name , unsigned int size) {
    if(polynomial.degree >= size) {
        print_pol(polynomial , name);
        return;
    }

    int blank , i=0;

    blank = size - polynomial.degree; 

    while(i<ceil((double)blank/2)) {
        if(!i) printf("   ");
        else printf("         ");
        i++;
    }
    print_pol(polynomial , name);
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
void add_pol(pol *res, pol A, pol B) {

    if(res->degree != MAX(A.degree , B.degree))
        change_degre_pol(res, MAX(A.degree, B.degree));

    for(int i=0 ; i<=res->degree ; i++) {
        if(i <= A.degree) {
            if(i <= B.degree) {
               res->coeffs[i] = A.coeffs[i] + B.coeffs[i];
            }
            else
                set_coeff_pol(*res, i, A.coeffs[i]);
        }
        else
            set_coeff_pol(*res, i, B.coeffs[i]);
    }
    while(!res->coeffs[res->degree] && res->degree>0) res->degree--;
}

/* ********************************************************************************************************************** */

/**
 * Subtract two polynomials
 * @param res
 * @param A
 * @param B
 */
void sub_pol(pol *res, pol A, pol B) {

    if(res->degree != MAX(A.degree , B.degree))
        change_degre_pol(res, MAX(A.degree, B.degree));

    for(int i=0 ; i<=res->degree ; i++) {
        if(i <= A.degree) {
            if(i <= B.degree) {
                res->coeffs[i] = A.coeffs[i] - B.coeffs[i];
            }
            else
                set_coeff_pol(*res, i, A.coeffs[i]);
        }
        else
            *(res->coeffs+i) = -1 * B.coeffs[i];
    }

    while(!res->coeffs[res->degree] && res->degree>0) res->degree--;
}

/* ********************************************************************************************************************** */

/**
 * Naive multiplication of two polynomials
 * @param res
 * @param A
 * @param B
 */
void mult_pol(pol *res, pol A, pol B) {
    if(is_zero_pol(A) || is_zero_pol(B)) {
        change_degre_pol(res , 0);
        set_coeff_pol(*res , 0 , 0);
        return;
    }
    change_degre_pol(res, A.degree + B.degree);
    long temp;

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

    while(!res->coeffs[res->degree] && res->degree>0) res->degree--;
}

/* ********************************************************************************************************************** */

void scalar_mult_pol(pol *res , pol A , long lambda) {
    change_degre_pol(res , A.degree);
    for(unsigned int i=0 ; i<res->degree ; i++)
        set_coeff_pol(*res , i , lambda*A.coeffs[i]);
}

/* ********************************************************************************************************************** */


int is_zero_pol(pol polynomial) {
    for(unsigned int i=0 ; i<=polynomial.degree ; i++)
        if (polynomial.coeffs[i]) return 0;
    return 1;
}

/* ********************************************************************************************************************** */


void euclide_div_pol(pol *Q , pol *R , pol A , pol B) {

   if(A.degree < B.degree) {
        perror("Euclidean division : A.degree < B.degree !!");
        return;
    }

    if(is_zero_pol(B)) {
        perror("Euclidean division : B can't be 0 !!");
        return;
    }   

    /* use double precision floating numbers to avoid round numbers */

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

    change_degre_pol(Q , sizeQ);
    change_degre_pol(R , sizeR);
    
    for(int i=0 ; i<=MAX(sizeR,sizeQ) ; i++) {
        if(i<=sizeR) R->coeffs[i] = (long)coeffsR[i];
        if(i<=sizeQ) Q->coeffs[i] = (long)coeffsQ[i];
    }

    free(coeffsA);
    free(coeffsB);
    free(coeffsT);
    free(coeffsR);
    free(coeffsQ);
}


/* ********************************************************************************************************************** */


void horner_eval(long *res , pol f , long x) {
    *res = f.coeffs[f.degree];
    for(unsigned int i=1 ; i<=f.degree ; i++) {
        *res *= x;
        *res += f.coeffs[f.degree - i];
    }
}

void horner_eval_multi(long *res , pol f , long *x , unsigned int nb_x) {
    for(unsigned int i=0 ; i<nb_x ; i++) horner_eval(res + i , f , *(x + i));
}


/* ********************************************************************************************************************** */

void derivate_pol(pol *derivate , pol A) {
    if(A.degree == 0) {
        init_pol(derivate , 0);
        return;
    }
    
    init_pol(derivate , A.degree-1);

    for(int i = 0 ; i <= derivate->degree ; i++) {
        derivate->coeffs[i] = A.coeffs[i+1] * (i+1);
    }
}


/* ********************************************************************************************************************** */

/*
* Transform the polynomial with the smaller coefficents possible considering the field Z/pZ
*/
void reduce_pol_ff(pol *A , long P) {
    for(int i = 0 ; i <= A->degree ; i++) {
        A->coeffs[i] %= P;
        A->coeffs[i] += P*(A->coeffs[i]<0);
    }
    for(;A->coeffs[A->degree] == 0 && A->degree > 0 ; --A->degree);
}

/* ********************************************************************************************************************** */

void euclide_div_pol_ff(pol *Q , pol *R , pol A , pol B , long P) {

     if(A.degree < B.degree) {
        perror("Euclidean division : A.degree < B.degree !!");
        return;
    }

    if(is_zero_pol(B)) {
        perror("Euclidean division : B can't be 0 !!");
        return;
    }  

    if(B.degree == 0) {
        init_pol(Q , A.degree);
        long alpha = modular_inverse(B.coeffs[0] , P);
        for(int i = 0 ; i <= Q->degree ; i++) {
            Q->coeffs[i] = A.coeffs[i]*alpha;
        }
        init_pol(R,0);
        return;
    }

    reduce_pol_ff(&A,P);
    reduce_pol_ff(&B,P);

    init_pol(R , A.degree);
    copy_pol(R,A);

    init_pol(Q , A.degree-B.degree);
    

    pol temp, R_temp; // temp = -(a/b)*x^(degR-degB)*B ; temp_R = R - temp
    init_pol(&temp,R->degree);
    init_pol(&R_temp , R->degree);

    long a,b,b_inv,t;
    b = B.coeffs[B.degree];
    b_inv = modular_inverse(b , P);    

    while(R->degree >= B.degree) {
        a = R->coeffs[R->degree];
        t = (a * b_inv) % P;

        Q->coeffs[R->degree-B.degree] += t; // head coefficient

        change_degre_pol(&temp,R->degree);
        for(int i = 0 ; i <=temp.degree ; i++) {
            if(i<R->degree-B.degree) temp.coeffs[i] = 0;
            else temp.coeffs[i] = -1 * t * B.coeffs[i-(R->degree-B.degree)];
        }
        add_pol(&R_temp,*R ,temp);

        reduce_pol_ff(&R_temp , P);
        
        copy_pol(R , R_temp);

        for(;R->coeffs[R->degree] == 0 && R->degree > 0 ; --R->degree);
    }

    for(;Q->coeffs[Q->degree] == 0 && Q->degree>0; --Q->degree);

    reduce_pol_ff(Q,P);
    reduce_pol_ff(R,P);  

    destroy_pol(temp);
    destroy_pol(R_temp);
}

/* ********************************************************************************************************************** */


void gcd_pol_ff(pol *gcd , pol A , pol B, long P) {

    // AA and BB are A and B such that degAA >= degBB
    pol AA, BB;

    if(A.degree < B.degree) {
       init_pol(&AA , B.degree);
       init_pol(&BB , A.degree);
       copy_pol(&AA , B);
       copy_pol(&BB , A);
    }
    else {
        init_pol(&AA , A.degree);
        init_pol(&BB , B.degree);
        copy_pol(&AA , A);
        copy_pol(&BB , B);
    }

    init_pol(gcd , MIN(AA.degree,BB.degree));

    pol R[3];

    init_pol(R , AA.degree);
    init_pol(R+1 , BB.degree);
    init_pol(R+2 , 0);

    reduce_pol_ff(&AA,P);
    reduce_pol_ff(&BB,P);

    copy_pol(R , AA);
    copy_pol(R+1 , BB);

    pol Q;
    init_pol(&Q , AA.degree-BB.degree);

    while(!is_zero_pol(R[1])) {
        euclide_div_pol_ff(&Q , R+2 , R[0] , R[1], P);
        copy_pol(R, R[1]);
        copy_pol(R+1 , R[2]);
        
    }

    copy_pol(gcd , R[0]);

}