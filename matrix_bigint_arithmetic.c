//
// Created by root on 30/10/18.
//

#include "matrix_bigint_arithmetic.h"

/**
 * Allocate the memory for the matrix values.
 * All the values are inited with 0
 * @param matrix
 * @param nb_line_p
 * @param nb_col_p
 */
void init_matrix_bigint(matrix_bigint *matrix , unsigned int nb_line_p , unsigned nb_col_p) {
   if(nb_line_p <= 0 || nb_col_p <= 0) {
       perror("Try to create a matrix without line or column");
       return;
   }

   matrix->values = (mpz_t *)calloc(sizeof(mpz_t) , nb_line_p*nb_col_p);
   for(int i=0 ; i<nb_line_p ; i++) {
       for(int j=0 ; j<nb_col_p ; j++) {
           mpz_init_set_d(matrix->values[i*nb_col_p + j],0);
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
void set_coeff_matrix_bigint(matrix_bigint matrix, unsigned int index_line , unsigned int index_col , mpz_t newvalue) {
    if(index_line < 0 || index_line >= matrix.nb_line || index_col < 0 || index_col >= matrix.nb_col) {
        perror("Matrix index not correct");
        return;
    }
    mpz_set(matrix.values[index_line*matrix.nb_col + index_col] , newvalue);

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
void set_coeff_matrix_bigint_d(matrix_bigint matrix, unsigned int index_line , unsigned int index_col , int newvalue) {
    if(index_line < 0 || index_line >= matrix.nb_line || index_col < 0 || index_col >= matrix.nb_col) {
        perror("Matrix index not correct");
        return;
    }
    mpz_set_d(matrix.values[index_line*matrix.nb_col + index_col] , newvalue);
}

/* ********************************************************************************************************************** */

/**
 * Set random values in the matrix between 0 and max
 * @param matrix
 * @param max
 */
void set_all_coeffs_random_matrix_bigint(matrix_bigint matrix , unsigned int max) {
    for(int i=0 ; i<matrix.nb_line ; i++) {
        for(int j=0 ; j<matrix.nb_col ; j++) {
            int val = (rand() % (2*max)) - max;
            mpz_set_d(matrix.values[i*matrix.nb_col + j] , val);
        }
    }
}

/* ********************************************************************************************************************** */

void destroy_matrix_bigint(matrix_bigint matrix) {
    for(int i=0 ; i<matrix.nb_line ; i++) {
        for(int j=0 ; j<matrix.nb_col ; j++) {
            mpz_clear(matrix.values[i*matrix.nb_col + j]);
        }
    }
    free(matrix.values);
}

/* ********************************************************************************************************************** */

void print_matrix_bigint(matrix_bigint matrix) {
    for(int i=0 ; i<matrix.nb_line ; i++) {
        for(int j=0 ; j<matrix.nb_col ; j++) {
            if(j == 0) printf("| ");

            if(mpz_cmp_d(matrix.values[i*matrix.nb_col + j],0) == 0) printf(" 0 ");
            if(mpz_cmp_d(matrix.values[i*matrix.nb_col + j],0) > 0) printf("+");
            if(mpz_cmp_d(matrix.values[i*matrix.nb_col + j],0) != 0) gmp_printf("%Zd ",matrix.values[i*matrix.nb_col + j]);
            if(mpz_cmp_d(matrix.values[i*matrix.nb_col + j],10) < 0  && mpz_cmp_d(matrix.values[i*matrix.nb_col + j],-10) > 0) printf(" ");
            if(j == matrix.nb_col - 1) printf("|\n");
        }
    }

    printf("\n");
}

/* ********************************************************************************************************************** */

void change_nb_line_matrix_bigint(matrix_bigint *matrix , unsigned int new_nb_line) {
    if(new_nb_line == matrix->nb_line) return;

    mpz_t *new_values = (mpz_t*)calloc(sizeof(mpz_t) , new_nb_line*matrix->nb_col);

    if(new_nb_line < matrix->nb_line) {
        for(int i=0 ; i<matrix->nb_line ; i++) {
            for(int j=0 ; j<matrix->nb_col ; j++) {
                if(i < new_nb_line)
                    mpz_init_set(new_values[i*matrix->nb_col + j] , matrix->values[i*matrix->nb_col + j]);

                mpz_clear(matrix->values[i*matrix->nb_col + j]);
            }
        }
    }
    else {
        for(int i=0 ; i<new_nb_line ; i++) {
            for(int j=0 ; j<matrix->nb_col ; j++) {
                if(i <matrix->nb_line) {
                    mpz_init_set(new_values[i*matrix->nb_col + j] , matrix->values[i*matrix->nb_col + j]);
                    mpz_clear(matrix->values[i*matrix->nb_col + j]);
                }
                else {
                    mpz_init_set_d(new_values[i*matrix->nb_col + j],0);
                }
            }
        }
    }

    free(matrix->values);
    matrix->values = new_values;
    matrix->nb_line = new_nb_line;
}

/* ********************************************************************************************************************** */

void change_nb_col_matrix_bigint(matrix_bigint *matrix , unsigned int new_nb_col) {
    if(new_nb_col == matrix->nb_col) return;

    mpz_t *newvalues = (mpz_t *)calloc(sizeof(mpz_t) , matrix->nb_line*new_nb_col);

    if(new_nb_col < matrix->nb_col) {
        for(int i=0 ; i<matrix->nb_line ; i++) {
            for(int j=0 ; j<matrix->nb_col ; j++) {
                if(j < new_nb_col)
                    mpz_init_set(newvalues[i*new_nb_col + j] , matrix->values[i*matrix->nb_col + j]);

                mpz_clear(matrix->values[i*matrix->nb_col + j]);
            }
        }
    }
    else {
        for(int i=0 ; i<matrix->nb_line ; i++) {
            for(int j=0 ; j<new_nb_col ; j++) {
                if(j < matrix->nb_col) {
                    mpz_init_set(newvalues[i*new_nb_col + j] , matrix->values[i*matrix->nb_col + j]);
                    mpz_clear(matrix->values[i*matrix->nb_col + j]);
                }
                else
                    mpz_init_set_d(newvalues[i*new_nb_col + j] , 0);
            }
        }
    }

    free(matrix->values);
    matrix->values = newvalues;
    matrix->nb_col = new_nb_col;
}

/* ********************************************************************************************************************** */

void change_dim_matrix_bigint(matrix_bigint *matrix , unsigned int new_nb_line , unsigned int new_nb_col) {
    if(new_nb_line != matrix->nb_line) change_nb_line_matrix_bigint(matrix , new_nb_line);
    if(new_nb_col != matrix->nb_col)   change_nb_col_matrix_bigint(matrix , new_nb_col);
}
