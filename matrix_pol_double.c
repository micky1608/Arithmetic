//
// Created by micky on 27/02/19.
//

#include "matrix_pol_double.h"

/**
 * Memory allocation for the polynomial matrix. I suppose that the memory for the polynomial has already been allocated.
 * It means that the function init_pol has been called for each polynomial before creating the matrix.
 * @param matrix_pol
 * @param nb_line
 * @param nb_col
 */
void init_matrix_pol_double(matrix_pol_double *matrix_pol , unsigned int nb_line , unsigned int nb_col) {
    if(nb_line <= 0 || nb_col<= 0) {
        perror("Try to create a matrix without line or column");
        return;
    }

    matrix_pol->nb_line = nb_line;
    matrix_pol->nb_col = nb_col;
    matrix_pol->values = (pol_double*)calloc(sizeof(pol_double) , nb_line*nb_col);
    for(int i=0 ; i<nb_line ; i++) {
        for(int j=0 ; j<nb_col ; j++) {
            init_pol_double(&MATRIX_P(matrix_pol,i,j) , 0);
        }
    }
}

/* ********************************************************************************************************************** */

/**
 * Replace one polynomial of the matrix with a new polynomial.
 * The memory allocated for the coefficients of the erased polynomial is freed.
*/
void setCoeff_matrix_pol_double(matrix_pol_double *matrix_pol , unsigned int line , unsigned int col , pol_double newvalue) {
    if(line < 0 || line >= matrix_pol->nb_line || col < 0 || col >= matrix_pol->nb_col) {
        perror("Matrix index not correct");
        return;
    }
    copy_pol_double(&MATRIX_P(matrix_pol,line,col) , newvalue);
}

/* ********************************************************************************************************************** */

void set_coeff_constant_matrix_pol_double(matrix_pol_double *matrix_pol , unsigned int line , unsigned int col , long value) {
    pol_double p;
    init_pol_double(&p , 0);
    set_coeff_pol_double(p , 0 , value);
    copy_pol_double(&MATRIX_P(matrix_pol,line,col) , p);
}

/* ********************************************************************************************************************** */

void set_allCoeff_matrix_pol_double(matrix_pol_double *matrix_pol , pol_double newvalue) {
     for(int i=0 ; i<matrix_pol->nb_line ; i++) {
        for(int j=0 ; j<matrix_pol->nb_col ; j++) {
            setCoeff_matrix_pol_double(matrix_pol , i , j , newvalue);
        }
     }
}

/* ********************************************************************************************************************** */

/**
 * Free the memory for the coefficients of each element of the matrix.
 * Then free the memory of each polynomial stucuture of the structure.
 */
void destroy_matrix_pol_double(matrix_pol_double matrix_pol) {
     for(int i=0 ; i<matrix_pol.nb_line ; i++) {
        for(int j=0 ; j<matrix_pol.nb_col ; j++) {
            destroy_pol_double(MATRIX(matrix_pol,i,j));
        }
     }
     free(matrix_pol.values);
}

/* ********************************************************************************************************************** */

void print_matrix_pol_double(matrix_pol_double matrix_pol , char *name) {
    int max_degree = 0;

    for(int j=0 ; j<matrix_pol.nb_col ; j++) {
        for(int i=0 ; i<matrix_pol.nb_line ; i++) {
            if(MATRIX(matrix_pol,i,j).degree > max_degree) {
                max_degree = MATRIX(matrix_pol,i,j).degree;
            }
        }
    }

    printf("%s :\n",name);
    for(int i=0 ; i<matrix_pol.nb_line ; i++) {
        for(int j=0 ; j<matrix_pol.nb_col ; j++) {
            if(j == 0) printf("| ");
            print_pol_double_center(MATRIX(matrix_pol,i,j) , NULL , max_degree);
            if(j == matrix_pol.nb_col - 1) printf(" |\n");
            else printf("     ");
        }
    }

    printf("\n");
}

/* ********************************************************************************************************************** */

void change_nb_line_matrix_pol_double(matrix_pol_double *matrix_pol , unsigned int new_nb_line) {
     if(new_nb_line == matrix_pol->nb_line) return;

    pol_double *new_values = (pol_double*)calloc(sizeof(pol_double) , new_nb_line*matrix_pol->nb_col);
    for(int i=0 ; i<new_nb_line*matrix_pol->nb_col ; i++) init_pol_double(&new_values[i] , 0);

    if(new_nb_line < matrix_pol->nb_line) {
        for(int i=0 ; i<matrix_pol->nb_line ; i++) {
            for(int j=0 ; j<matrix_pol->nb_col ; j++) {
                if(i < new_nb_line) {
                    copy_pol_double(&new_values[i*matrix_pol->nb_col + j] , MATRIX_P(matrix_pol,i,j));
                }
            }
        }
    }
    else {
        for(int i=0 ; i<new_nb_line ; i++) {
            for(int j=0 ; j<matrix_pol->nb_col ; j++) {
                if(i <matrix_pol->nb_line) {
                    copy_pol_double(&new_values[i*matrix_pol->nb_col + j] , MATRIX_P(matrix_pol,i,j));
                }
                else {
                    change_degre_pol_double(&new_values[i*matrix_pol->nb_col + j] , 0);
                    set_coeff_pol_double(new_values[i*matrix_pol->nb_col + j] , 0 , 0);
                }
            }
        }
    }

    for(int i=0 ; i<matrix_pol->nb_line ; i++) {
        for(int j=0 ; j<matrix_pol->nb_col ; j++) {
            destroy_pol_double(MATRIX_P(matrix_pol,i,j));
        }
     }
    free(matrix_pol->values);
    matrix_pol->values = new_values;
    matrix_pol->nb_line = new_nb_line;
}

/* ********************************************************************************************************************** */

void change_nb_col_matrix_pol_double(matrix_pol_double *matrix_pol , unsigned int new_nb_col) {
     if(new_nb_col == matrix_pol->nb_col) return;

    pol_double *newvalues = (pol_double*)calloc(sizeof(pol_double) , matrix_pol->nb_line*new_nb_col);
    for(int i=0 ; i<matrix_pol->nb_line*new_nb_col ; i++) init_pol_double(&newvalues[i] , 0);

    if(new_nb_col < matrix_pol->nb_col) {
        for(int i=0 ; i<matrix_pol->nb_line ; i++) {
            for(int j=0 ; j<matrix_pol->nb_col ; j++) {
                if(j < new_nb_col) {
                    //newvalues[i*new_nb_col + j] = matrix->values[i*matrix->nb_col + j];
                    copy_pol_double(&newvalues[i*new_nb_col + j] , MATRIX_P(matrix_pol,i,j));
                }
            }
        }
    }
    else {
        for(int i=0 ; i<matrix_pol->nb_line ; i++) {
            for(int j=0 ; j<new_nb_col ; j++) {
                if(j < matrix_pol->nb_col) {
                    //newvalues[i*new_nb_col + j] = matrix->values[i*matrix->nb_col + j];
                    copy_pol_double(&newvalues[i*new_nb_col + j] , MATRIX_P(matrix_pol,i,j));
                }
                else {
                    change_degre_pol_double(&newvalues[i*new_nb_col + j] , 0);
                    set_coeff_pol_double(newvalues[i*new_nb_col + j] , 0 , 0);
                }
            }
        }
    }

    for(int i=0 ; i<matrix_pol->nb_line ; i++) {
        for(int j=0 ; j<matrix_pol->nb_col ; j++) {
            destroy_pol_double(MATRIX_P(matrix_pol,i,j));
        }
     }
    free(matrix_pol->values);
    matrix_pol->values = newvalues;
    matrix_pol->nb_col = new_nb_col;
}

/* ********************************************************************************************************************** */

void change_dim_matrix_pol_double(matrix_pol_double *matrix_pol , unsigned int new_nb_line , unsigned int new_nb_col) {
    if(new_nb_line != matrix_pol->nb_line) change_nb_line_matrix_pol_double(matrix_pol , new_nb_line);
    if(new_nb_col != matrix_pol->nb_col)   change_nb_col_matrix_pol_double(matrix_pol , new_nb_col);
}

/* ********************************************************************************************************************** */

/**
 * Create a polynomial identity matrix.
 */
void identity_matrix_pol_double(matrix_pol_double *id, unsigned int size) {
    init_matrix_pol_double(id, size, size);

    for(int i=0 ; i<size ; i++) {
        for(int j=0 ; j<size ; j++) {
            if(i == j) set_coeff_pol_double(MATRIX_P(id,i,j) , 0 , 1);
        }
    }
}

/* ********************************************************************************************************************** */

void add_matrix_pol_double(matrix_pol_double *res , matrix_pol_double A , matrix_pol_double B) {
     if(A.nb_line != B.nb_line || A.nb_col != B.nb_col) {
        perror("Add matrix : error dimensions");
        return;
    }

    change_dim_matrix_pol_double(res , A.nb_line , A.nb_col);

    for(unsigned int i=0 ; i<res->nb_line ; i++) {
        for(unsigned int j=0 ; j<res->nb_col ; j++)
            add_pol_double(&MATRIX_P(res,i,j) , MATRIX(A,i,j) , MATRIX(B,i,j));
    }
}

/* ********************************************************************************************************************** */

void mul_matrix_pol_double(matrix_pol_double *res , matrix_pol_double A , matrix_pol_double B) {
     if(A.nb_col != B.nb_line) {
        perror("Error dimension mul_matrix_pol");
        exit(EXIT_FAILURE);
    }

    init_matrix_pol_double(res , A.nb_line , B.nb_col);

    pol_double temp_product , temp_add;
    init_pol_double(&temp_product , A.values[0].degree + B.values[0].degree);
    init_pol_double(&temp_add , 0);

    for(int i=0 ; i<res->nb_line ; i++) {
        for(int j=0 ; j<res->nb_col ; j++) {
            for(int k=0 ; k<A.nb_line ; k++) {
                mult_pol_double(&temp_product , MATRIX(A,i,k) , MATRIX(B,k,j));
                add_pol_double(&temp_add , temp_product , MATRIX_P(res,i,j));
                copy_pol_double(&MATRIX_P(res,i,j) , temp_add);
            }
        }
    }

    destroy_pol_double(temp_product);
    destroy_pol_double(temp_add);
}

/* ********************************************************************************************************************** */

void copy_matrix_pol_double(matrix_pol_double *DEST , matrix_pol_double SRC) {
    change_dim_matrix_pol_double(DEST , SRC.nb_line , SRC.nb_col);

    for(int i=0 ; i<SRC.nb_line ; i++) {
        for (int j = 0; j < SRC.nb_col; j++) {
            copy_pol_double(&MATRIX_P(DEST,i,j) , MATRIX(SRC,i,j));
        }
    }
}
