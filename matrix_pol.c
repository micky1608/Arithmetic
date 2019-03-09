//
// Created by micky on 27/02/19.
//

#include "matrix_pol.h"

/**
 * Memory allocation for the polynomial matrix. I suppose that the memory for the polynomial has already been allocated.
 * It means that the function init_pol has been called for each polynomial before creating the matrix.
 * @param matrix_pol
 * @param nb_line
 * @param nb_col
 */
void init_matrix_pol(matrix_pol *matrix_pol , unsigned int nb_line , unsigned int nb_col) {
    if(nb_line <= 0 || nb_col<= 0) {
        perror("Try to create a matrix without line or column");
        return;
    }

    matrix_pol->nb_line = nb_line;
    matrix_pol->nb_col = nb_col;
    matrix_pol->values = (pol*)calloc(sizeof(pol) , nb_line*nb_col);
    for(int i=0 ; i<nb_line ; i++) {
        for(int j=0 ; j<nb_col ; j++) {
            init_pol(&MATRIX_P(matrix_pol,i,j) , 0);
        }
    }
}

/* ********************************************************************************************************************** */

/**
 * Replace one polynomial of the matrix with a new polynomial.
 * The memory allocated for the coefficients of the erased polynomial is freed.
*/
void setCoeff_matrix_pol(matrix_pol *matrix_pol , unsigned int line , unsigned int col , pol newvalue) {
    if(line < 0 || line >= matrix_pol->nb_line || col < 0 || col >= matrix_pol->nb_col) {
        perror("Matrix index not correct");
        return;
    }
    copy_pol(&MATRIX_P(matrix_pol,line,col) , newvalue);
}

/* ********************************************************************************************************************** */

void set_allCoeff_matrix_pol(matrix_pol *matrix_pol , pol newvalue) {
     for(int i=0 ; i<matrix_pol->nb_line ; i++) {
        for(int j=0 ; j<matrix_pol->nb_col ; j++) {
            setCoeff_matrix_pol(matrix_pol , i , j , newvalue);
        }
     }
}

/* ********************************************************************************************************************** */

/**
 * Free the memory for the coefficients of each element of the matrix.
 * Then free the memory of each polynomial stucuture of the structure.
 */
void destroy_matrix_pol(matrix_pol matrix_pol) {
     for(int i=0 ; i<matrix_pol.nb_line ; i++) {
        for(int j=0 ; j<matrix_pol.nb_col ; j++) {
            destroy_pol(MATRIX(matrix_pol,i,j));
        }
     }
     free(matrix_pol.values);
}

/* ********************************************************************************************************************** */

void print_matrix_pol(matrix_pol matrix_pol , char *name) {
    int max_degree_col[matrix_pol.nb_col];
    memset(max_degree_col , 0 , matrix_pol.nb_col * sizeof(int));

    for(int j=0 ; j<matrix_pol.nb_col ; j++) {
        for(int i=0 ; i<matrix_pol.nb_line ; i++) {
            if(MATRIX(matrix_pol,i,j).degree > max_degree_col[j]) {
                max_degree_col[j] = MATRIX(matrix_pol,i,j).degree;
            }
        }
    }

    printf("%s :\n",name);
    for(int i=0 ; i<matrix_pol.nb_line ; i++) {
        for(int j=0 ; j<matrix_pol.nb_col ; j++) {
            if(j == 0) printf("| ");
            print_pol_center(MATRIX(matrix_pol,i,j) , NULL , max_degree_col[j]);
            if(j == matrix_pol.nb_col - 1) printf(" |\n");
            else printf("     ");
        }
    }

    printf("\n");
}

/* ********************************************************************************************************************** */

void change_nb_line_matrix_pol(matrix_pol *matrix_pol , unsigned int new_nb_line) {
     if(new_nb_line == matrix_pol->nb_line) return;

    pol *new_values = (pol*)calloc(sizeof(pol) , new_nb_line*matrix_pol->nb_col);
    for(int i=0 ; i<new_nb_line*matrix_pol->nb_col ; i++) init_pol(&new_values[i] , 0);

    if(new_nb_line < matrix_pol->nb_line) {
        for(int i=0 ; i<matrix_pol->nb_line ; i++) {
            for(int j=0 ; j<matrix_pol->nb_col ; j++) {
                if(i < new_nb_line) {
                    //new_values[i*matrix->nb_col + j] = MATRIX_P(matrix,i,j);
                    copy_pol(&new_values[i*matrix_pol->nb_col + j] , MATRIX_P(matrix_pol,i,j));
                }
            }
        }
    }
    else {
        for(int i=0 ; i<new_nb_line ; i++) {
            for(int j=0 ; j<matrix_pol->nb_col ; j++) {
                if(i <matrix_pol->nb_line) {
                    //new_values[i*matrix_pol->nb_col + j] = MATRIX_P(matrix,i,j);
                    copy_pol(&new_values[i*matrix_pol->nb_col + j] , MATRIX_P(matrix_pol,i,j));
                }
                else {
                    change_degre_pol(&new_values[i*matrix_pol->nb_col + j] , 0);
                    set_coeff_pol(new_values[i*matrix_pol->nb_col + j] , 0 , 0);
                }
            }
        }
    }

    for(int i=0 ; i<matrix_pol->nb_line ; i++) {
        for(int j=0 ; j<matrix_pol->nb_col ; j++) {
            destroy_pol(MATRIX_P(matrix_pol,i,j));
        }
     }
    free(matrix_pol->values);
    matrix_pol->values = new_values;
    matrix_pol->nb_line = new_nb_line;
}

/* ********************************************************************************************************************** */

void change_nb_col_matrix_pol(matrix_pol *matrix_pol , unsigned int new_nb_col) {
     if(new_nb_col == matrix_pol->nb_col) return;

    pol *newvalues = (pol*)calloc(sizeof(pol) , matrix_pol->nb_line*new_nb_col);
    for(int i=0 ; i<matrix_pol->nb_line*new_nb_col ; i++) init_pol(&newvalues[i] , 0);

    if(new_nb_col < matrix_pol->nb_col) {
        for(int i=0 ; i<matrix_pol->nb_line ; i++) {
            for(int j=0 ; j<matrix_pol->nb_col ; j++) {
                if(j < new_nb_col) {
                    //newvalues[i*new_nb_col + j] = matrix->values[i*matrix->nb_col + j];
                    copy_pol(&newvalues[i*new_nb_col + j] , MATRIX_P(matrix_pol,i,j));
                }
            }
        }
    }
    else {
        for(int i=0 ; i<matrix_pol->nb_line ; i++) {
            for(int j=0 ; j<new_nb_col ; j++) {
                if(j < matrix_pol->nb_col) {
                    //newvalues[i*new_nb_col + j] = matrix->values[i*matrix->nb_col + j];
                    copy_pol(&newvalues[i*new_nb_col + j] , MATRIX_P(matrix_pol,i,j));
                }
                else {
                    change_degre_pol(&newvalues[i*new_nb_col + j] , 0);
                    set_coeff_pol(newvalues[i*new_nb_col + j] , 0 , 0);
                }
            }
        }
    }

    for(int i=0 ; i<matrix_pol->nb_line ; i++) {
        for(int j=0 ; j<matrix_pol->nb_col ; j++) {
            destroy_pol(MATRIX_P(matrix_pol,i,j));
        }
     }
    free(matrix_pol->values);
    matrix_pol->values = newvalues;
    matrix_pol->nb_col = new_nb_col;
}

/* ********************************************************************************************************************** */

void change_dim_matrix_pol(matrix_pol *matrix_pol , unsigned int new_nb_line , unsigned int new_nb_col) {
    if(new_nb_line != matrix_pol->nb_line) change_nb_line_matrix_pol(matrix_pol , new_nb_line);
    if(new_nb_col != matrix_pol->nb_col)   change_nb_col_matrix_pol(matrix_pol , new_nb_col);
}

/* ********************************************************************************************************************** */

void mul_matrix_pol(matrix_pol *res , matrix_pol A , matrix_pol B);

/* ********************************************************************************************************************** */

void copy_matrix_pol(matrix_pol *DEST , matrix_pol SRC);
