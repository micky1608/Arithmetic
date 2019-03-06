//
// Created by micky on 27/02/19.
//

#include "matrix_pol.h"

/**
 * Memory allocation for the polynomial matrix. I suppose that the memory for the polynomial has already been allocated.
 * This function reserves memory to store the polynomials structures.
 * @param matrix_pol
 * @param nb_line
 * @param nb_col
 */
void init_matrix_pol(matrix_pol *matrix_pol , unsigned int nb_line , unsigned int nb_col) {
    if(nb_line <= 0 || nb_col<= 0) {
        perror("Try to create a matrix without line or column");
        return;
    }

    matrix_pol->values = (pol*)calloc(sizeof(pol) , nb_line*nb_col);
    for(int i=0 ; i<nb_line ; i++) {
        for(int j=0 ; j<nb_col ; j++) {
            set_all_coeffs_to_pol(matrix_pol->values[i*nb_col + j] , 0);
        }
    }
    matrix_pol->nb_line = nb_line;
    matrix_pol->nb_col = nb_col;
}

/* ********************************************************************************************************************** */

/**
 * Call the above function and set the coefficients of all the polynomials to 0
 * @param matrix_pol
 * @param nb_line
 * @param nb_col
 */
void init_matrix_zero_pol(matrix_pol *matrix_pol , unsigned int nb_line , unsigned int nb_col) {
    init_matrix_pol(matrix_pol , nb_line , nb_col);
    for(int i=0 ; i<nb_line ; i++) {
        for(int j=0 ; j<nb_col ; j++) {
            set_all_coeffs_to_pol(matrix_pol->values[i*nb_col + j] , 0);
        }
    }
}

/* ********************************************************************************************************************** */

void setCoeff_matrix_pol(matrix_pol *matrix_pol , unsigned int line , unsigned int col , pol newvalue) {
    if(line < 0 || line >= matrix_pol->nb_line || col < 0 || col >= matrix_pol->nb_col) {
        perror("Matrix index not correct");
        return;
    }
    free(MATRIX_P(matrix_pol,line,col).coeffs);
    MATRIX_P(matrix_pol,line,col) = newvalue;
}

/* ********************************************************************************************************************** */

void setCoeff_matrix_pol_array(matrix_pol *matrix_pol , pol *coeffs , unsigned int sizeArray) {
    if(matrix_pol->nb_line*matrix_pol->nb_col != sizeArray) {
        perror("SetCoeffArray : Error array size not correct");
        exit(EXIT_FAILURE);
    }

    for(unsigned int i=0 ; i<sizeArray ; i++)
        matrix_pol->values[i] = coeffs[i];
}

/* ********************************************************************************************************************** */

void setAllCoeff_matrix_pol(matrix_pol *matrix_pol , long value);

/* ********************************************************************************************************************** */

void destroy_matrix_pol(matrix_pol matrix_pol);

/* ********************************************************************************************************************** */

void print_matrix_pol(matrix_pol matrix_pol , char *name);

/* ********************************************************************************************************************** */

void change_nb_line_matrix_pol(matrix_pol *matrix_pol , unsigned int new_nb_line);

/* ********************************************************************************************************************** */

void change_nb_col_matrix_pol(matrix_pol *matrix_pol , unsigned int new_nb_col);

/* ********************************************************************************************************************** */

void change_dim_matrix_pol(matrix_pol *matrix_pol , unsigned int new_nb_line , unsigned int new_nb_col);

/* ********************************************************************************************************************** */

void mul_matrix_pol(matrix_pol *res , matrix_pol A , matrix_pol B);

/* ********************************************************************************************************************** */

void copy_matrix_pol(matrix_pol *DEST , matrix_pol SRC);
