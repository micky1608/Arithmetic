//
// Created by root on 30/10/18.
//

#include "matrix_bigQ_arithmetic.h"

/**
 * Allocate the memory for the matrix values.
 * All the values are inited with 0
 * @param matrix
 * @param nb_line_p
 * @param nb_col_p
 */
void init_matrix_bigQ(matrix_bigQ *matrix , unsigned int nb_line_p , unsigned nb_col_p) {
   if(nb_line_p <= 0 || nb_col_p <= 0) {
       perror("Try to create a matrix without line or column");
       return;
   }

   matrix->values = (mpq_t *)calloc(sizeof(mpq_t) , nb_line_p*nb_col_p);
   for(int i=0 ; i<nb_line_p ; i++) {
       for(int j=0 ; j<nb_col_p ; j++) {
           mpq_init(matrix->values[i*nb_col_p + j]);
           mpq_set_d(matrix->values[i*nb_col_p + j],0);
       }
   }

   matrix->nb_line = nb_line_p;
   matrix->nb_col = nb_col_p;
}

/* ********************************************************************************************************************** */

/**
 * Modify the value at the position (index_line,index_col) in the matrix.
 * WARNING : the numbering of lines and columns starts with 0
 * @param matrix
 * @param index_line
 * @param index_col
 * @param newvalue
 */
void set_coeff_matrix_bigQ(matrix_bigQ matrix, unsigned int index_line , unsigned int index_col , mpq_t newvalue) {
    if(index_line < 0 || index_line >= matrix.nb_line || index_col < 0 || index_col >= matrix.nb_col) {
        perror("Matrix index not correct");
        return;
    }
    mpq_set(matrix.values[index_line*matrix.nb_col + index_col] , newvalue);

}

/* ********************************************************************************************************************** */

/**
* Modify the value at the position (index_line,index_col) in the matrix.
* WARNING : the numbering of lines and columns starts with 0
* @param matrix
* @param index_line
* @param index_col
* @param newvalue
*/
void set_coeff_matrix_bigQ_d(matrix_bigQ matrix, unsigned int index_line , unsigned int index_col , int numerator , unsigned int denominator) {
    if(index_line < 0 || index_line >= matrix.nb_line || index_col < 0 || index_col >= matrix.nb_col) {
        perror("Matrix index not correct");
        return;
    }

    if(!denominator) {
        perror("Denominator can't be 0");
        return;
    }

    mpq_set_si(matrix.values[index_line*matrix.nb_col + index_col] , numerator, denominator);
}

/* ********************************************************************************************************************** */

/**
 * Set random values in the matrix between 0 and max
 * @param matrix
 * @param max
 */
void set_all_coeffs_random_matrix_bigQ(matrix_bigQ matrix , unsigned int max) {
    for(int i=0 ; i<matrix.nb_line ; i++) {
        for(int j=0 ; j<matrix.nb_col ; j++) {
            int numerator = (rand() % (2*max)) - max;
            unsigned int denominator;
            do {
                denominator = rand() % max;
            } while(!denominator);

            mpq_set_si(matrix.values[i*matrix.nb_col + j] , numerator , denominator);
        }
    }
}

/* ********************************************************************************************************************** */

void destroy_matrix_bigQ(matrix_bigQ matrix) {
    for(int i=0 ; i<matrix.nb_line ; i++) {
        for(int j=0 ; j<matrix.nb_col ; j++) {
            mpq_clear(matrix.values[i*matrix.nb_col + j]);
        }
    }
    free(matrix.values);
}

/* ********************************************************************************************************************** */

void print_matrix_bigQ(matrix_bigQ matrix) {
    for(int i=0 ; i<matrix.nb_line ; i++) {
        for(int j=0 ; j<matrix.nb_col ; j++) {
            if(j == 0) printf("| ");

            if(mpq_cmp_si(matrix.values[i*matrix.nb_col + j],0,1) == 0) printf(" 0 ");
            if(mpq_cmp_si(matrix.values[i*matrix.nb_col + j],0,1) > 0) printf("+");
            if(mpq_cmp_si(matrix.values[i*matrix.nb_col + j],0,1) != 0) gmp_printf("%Qx ",matrix.values[i*matrix.nb_col + j]);
            if(mpq_cmp_si(matrix.values[i*matrix.nb_col + j],10,1) < 0  && mpq_cmp_si(matrix.values[i*matrix.nb_col + j],-10,1) > 0) printf(" ");
            if(j == matrix.nb_col - 1) printf("|\n");
        }
    }

    printf("\n");
}

/* ********************************************************************************************************************** */

void change_nb_line_matrix_bigQ(matrix_bigQ *matrix , unsigned int new_nb_line) {
    if(new_nb_line == matrix->nb_line) return;

    mpq_t *new_values = (mpq_t*)calloc(sizeof(mpq_t) , new_nb_line*matrix->nb_col);

    if(new_nb_line < matrix->nb_line) {
        for(int i=0 ; i<matrix->nb_line ; i++) {
            for(int j=0 ; j<matrix->nb_col ; j++) {
                if(i < new_nb_line) {
                    mpq_init(new_values[i*matrix->nb_col + j]);
                    mpq_set(new_values[i*matrix->nb_col + j] , matrix->values[i*matrix->nb_col + j]);
                }

                mpq_clear(matrix->values[i*matrix->nb_col + j]);
            }
        }
    }
    else {
        for(int i=0 ; i<new_nb_line ; i++) {
            for(int j=0 ; j<matrix->nb_col ; j++) {
                if(i <matrix->nb_line) {
                    mpq_init(new_values[i*matrix->nb_col + j]);
                    mpq_set(new_values[i*matrix->nb_col + j] , matrix->values[i*matrix->nb_col + j]);
                    mpq_clear(matrix->values[i*matrix->nb_col + j]);
                }
                else {
                    mpq_init(new_values[i*matrix->nb_col + j]);
                    mpq_set_d(new_values[i*matrix->nb_col + j],0);
                }
            }
        }
    }

    free(matrix->values);
    matrix->values = new_values;
    matrix->nb_line = new_nb_line;
}

/* ********************************************************************************************************************** */

void change_nb_col_matrix_bigQ(matrix_bigQ *matrix , unsigned int new_nb_col) {
    if(new_nb_col == matrix->nb_col) return;

    mpq_t *newvalues = (mpq_t *)calloc(sizeof(mpq_t) , matrix->nb_line*new_nb_col);

    if(new_nb_col < matrix->nb_col) {
        for(int i=0 ; i<matrix->nb_line ; i++) {
            for(int j=0 ; j<matrix->nb_col ; j++) {
                if(j < new_nb_col) {
                    mpq_init(newvalues[i*new_nb_col + j]);
                    mpq_set(newvalues[i*new_nb_col + j] , matrix->values[i*matrix->nb_col + j]);
                }

                mpq_clear(matrix->values[i*matrix->nb_col + j]);
            }
        }
    }
    else {
        for(int i=0 ; i<matrix->nb_line ; i++) {
            for(int j=0 ; j<new_nb_col ; j++) {
                if(j < matrix->nb_col) {
                    mpq_init(newvalues[i*new_nb_col + j]);
                    mpq_set(newvalues[i*new_nb_col + j] , matrix->values[i*matrix->nb_col + j]);
                    mpq_clear(matrix->values[i*matrix->nb_col + j]);
                }
                else {
                    mpq_init(newvalues[i*new_nb_col + j]);
                    mpq_set_d(newvalues[i*new_nb_col + j] , 0);
                }

            }
        }
    }

    free(matrix->values);
    matrix->values = newvalues;
    matrix->nb_col = new_nb_col;
}

/* ********************************************************************************************************************** */

void change_dim_matrix_bigQ(matrix_bigQ *matrix , unsigned int new_nb_line , unsigned int new_nb_col) {
    if(new_nb_line != matrix->nb_line) change_nb_line_matrix_bigQ(matrix , new_nb_line);
    if(new_nb_col != matrix->nb_col)   change_nb_col_matrix_bigQ(matrix , new_nb_col);
}

/* ********************************************************************************************************************** */

void add_matrix_bigQ(matrix_bigQ *res , matrix_bigQ A , matrix_bigQ B) {
    if(A.nb_line != B.nb_line || A.nb_col != B.nb_col) {
        perror("Try to add two matrix with different dimensions");
        return;
    }

    change_dim_matrix_bigQ(res , A.nb_line , A.nb_col);

    for(int i=0 ; i<res->nb_line ; i++) {
        for(int j=0 ; j<res->nb_col ; j++) {
            mpq_add(res->values[i*res->nb_col + j] , A.values[i*res->nb_col + j] , B.values[i*res->nb_col + j]);
        }
    }
}

/* ********************************************************************************************************************** */

void sub_matrix_bigQ(matrix_bigQ *res , matrix_bigQ A , matrix_bigQ B) {

    matrix_bigQ B_neg;
    init_matrix_bigQ(&B_neg , B.nb_line , B.nb_col);
    scalar_mult_matrix_bigQ(&B_neg , B , -1);
    add_matrix_bigQ(res , A , B_neg);
    destroy_matrix_bigQ(B_neg);
}

/* ********************************************************************************************************************** */

void scalar_mult_matrix_bigQ(matrix_bigQ *res , matrix_bigQ A , long lambda) {
    change_dim_matrix_bigQ(res , A.nb_line , A.nb_col);

    mpz_t new_numerator;
    mpz_init(new_numerator);

    for(int i=0 ; i<res->nb_line ; i++) {
        for (int j = 0; j < res->nb_col; j++) {
            mpq_get_num(new_numerator, A.values[i*res->nb_col + j]);
            mpz_mul_si(new_numerator , new_numerator , lambda);
            mpq_set_num(res->values[i*res->nb_col + j] , new_numerator);
        }
    }

    mpz_clear(new_numerator);
}

/* ********************************************************************************************************************** */

void mult_matrix_bigQ(matrix_bigQ *res , matrix_bigQ A , matrix_bigQ B) {
    if(A.nb_col != B.nb_line) {
        perror("Try to multiply matrix with incorrect dimensions");
        return;
    }

    mpq_t temp;
    mpq_init(temp);
    change_dim_matrix_bigQ(res , A.nb_line , B.nb_col);

    for(int i=0 ; i<res->nb_line ; i++) {
        for(int j=0 ; j<res->nb_col ; j++) {
            mpq_set_d(res->values[i*res->nb_col + j] , 0);

            for(int k=0 ; k<A.nb_col ; k++) {
                mpq_mul(temp , A.values[i*A.nb_col + k] , B.values[k*B.nb_col + j]);
                mpq_add(res->values[i*res->nb_col + j] , res->values[i*res->nb_col + j] , temp);
            }
        }
    }

    mpq_clear(temp);
}
