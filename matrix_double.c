//
// Created by root on 07/12/18.
//

#include "matrix_double.h"

void init_matrix_double(matrix_double *matrixDouble , unsigned int nb_line_p , unsigned int nb_col_p) {
    if(nb_line_p <= 0 || nb_col_p <= 0) {
        perror("Try to create a matrix without line or column");
        return;
    }

    matrixDouble->values = (double*)calloc(sizeof(double) , nb_line_p*nb_col_p);
    for(int i=0 ; i<nb_line_p ; i++) {
        for(int j=0 ; j<nb_col_p ; j++) {
            matrixDouble->values[i*nb_col_p + j] = 0;
        }
    }

    matrixDouble->nb_line = nb_line_p;
    matrixDouble->nb_col = nb_col_p;
}

/* ********************************************************************************************************************** */


void setCoeff_matrix_double(matrix_double *matrixDouble , unsigned int line , unsigned int col , double newvalue) {
    if(line < 0 || line >= matrixDouble->nb_line || col < 0 || col >= matrixDouble->nb_col) {
        perror("Matrix index not correct");
        return;
    }
    matrixDouble->values[line*matrixDouble->nb_col + col] = newvalue;
}

/* ********************************************************************************************************************** */

void destroy_matrix_double(matrix_double matrixDouble) {
    free(matrixDouble.values);
}

/* ********************************************************************************************************************** */

void print_matrix_double(matrix_double matrixDouble , char *name) {
    printf("%s :\n",name);
    for(int i=0 ; i<matrixDouble.nb_line ; i++) {
        for(int j=0 ; j<matrixDouble.nb_col ; j++) {
            if(j == 0) printf("| ");

            if(matrixDouble.values[i*matrixDouble.nb_col + j] == 0) printf(" 0.00 ");
            if(matrixDouble.values[i*matrixDouble.nb_col + j] > 0) printf("+");
            if(matrixDouble.values[i*matrixDouble.nb_col + j] != 0) printf("%.2lf ",matrixDouble.values[i*matrixDouble.nb_col + j]);
            if(matrixDouble.values[i*matrixDouble.nb_col + j] < 10  && matrixDouble.values[i*matrixDouble.nb_col + j] > -10) printf(" ");
            if(j == matrixDouble.nb_col - 1) printf("|\n");
        }
    }

    printf("\n");
}

/* ********************************************************************************************************************** */

void change_nb_line_matrix_double(matrix_double *matrixDouble , unsigned int new_nb_line) {
    if(new_nb_line == matrixDouble->nb_line) return;

    double *new_values = (double*)calloc(sizeof(double) , new_nb_line*matrixDouble->nb_col);

    if(new_nb_line < matrixDouble->nb_line) {
        for(int i=0 ; i<matrixDouble->nb_line ; i++) {
            for(int j=0 ; j<matrixDouble->nb_col ; j++) {
                if(i < new_nb_line) {
                    new_values[i*matrixDouble->nb_col + j] = matrixDouble->values[i*matrixDouble->nb_col + j];
                }
            }
        }
    }
    else {
        for(int i=0 ; i<new_nb_line ; i++) {
            for(int j=0 ; j<matrixDouble->nb_col ; j++) {
                if(i <matrixDouble->nb_line) {
                    new_values[i*matrixDouble->nb_col + j] = matrixDouble->values[i*matrixDouble->nb_col + j];
                }
                else {
                    new_values[i*matrixDouble->nb_col + j] = 0;
                }
            }
        }
    }
    free(matrixDouble->values);
    matrixDouble->values = new_values;
    matrixDouble->nb_line = new_nb_line;
}

/* ********************************************************************************************************************** */

void change_nb_col_matrix_double(matrix_double *matrixDouble , unsigned int new_nb_col) {
    if(new_nb_col == matrixDouble->nb_col) return;

    double *newvalues = (double*)calloc(sizeof(double) , matrixDouble->nb_line*new_nb_col);

    if(new_nb_col < matrixDouble->nb_col) {
        for(int i=0 ; i<matrixDouble->nb_line ; i++) {
            for(int j=0 ; j<matrixDouble->nb_col ; j++) {
                if(j < new_nb_col) {
                    newvalues[i*new_nb_col + j] = matrixDouble->values[i*matrixDouble->nb_col + j];
                }
            }
        }
    }
    else {
        for(int i=0 ; i<matrixDouble->nb_line ; i++) {
            for(int j=0 ; j<new_nb_col ; j++) {
                if(j < matrixDouble->nb_col) {
                    newvalues[i*new_nb_col + j] = matrixDouble->values[i*matrixDouble->nb_col + j];
                }
                else {
                    newvalues[i*new_nb_col + j] = 0;
                }
            }
        }
    }

    free(matrixDouble->values);
    matrixDouble->values = newvalues;
    matrixDouble->nb_col = new_nb_col;
}

/* ********************************************************************************************************************** */

void change_dim_matrix_double(matrix_double *matrixDouble , unsigned int new_nb_line , unsigned int new_nb_col) {
    if(new_nb_line != matrixDouble->nb_line) change_nb_line_matrix_double(matrixDouble , new_nb_line);
    if(new_nb_col != matrixDouble->nb_col)   change_nb_col_matrix_double(matrixDouble , new_nb_col);
}

/* ********************************************************************************************************************** */

void copy_matrix_double(matrix_double *DEST , matrix_double SRC) {
    change_dim_matrix_double(DEST , SRC.nb_line , SRC.nb_col);

    for(int i=0 ; i<SRC.nb_line ; i++) {
        for (int j = 0; j < SRC.nb_col; j++) {
            DEST->values[i*SRC.nb_col + j] = SRC.values[i*SRC.nb_col + j];
        }
    }
}

/* ********************************************************************************************************************** */

void identity_matrix_double(matrix_double *id, unsigned int size) {
    change_dim_matrix_double(id, size, size);

    for(int i=0 ; i<size ; i++) {
        for(int j=0 ; j<size ; j++) {
            if(i == j)
                id->values[i*size +j] = 1;
            else
                id->values[i*size +j] = 0;
        }
    }
}

/* ********************************************************************************************************************** */

void LU_decomposition_matrix_double(matrix_double *L , matrix_double *U , matrix_double A) {

    if(L != NULL)
        identity_matrix_double(L , A.nb_line); // L = Id
    copy_matrix_double(U , A);             // U = A

    double coeff, temp;

    for(int j=0 ; j<A.nb_col ; j++) {
        for(int i=j+1 ; i<A.nb_line ; i++) {

            coeff = -1 * U->values[i*U->nb_col + j] / U->values[j*U->nb_col + j];

            for(int k=0 ; k<A.nb_col ; k++) {
                temp = U->values[j*U->nb_col+k] * coeff;

                if(L != NULL)
                    L->values[i*L->nb_col+j] = coeff;

                U->values[i*U->nb_col+k] += temp;
            }
        }
    }
}

/* ********************************************************************************************************************** */

void swap_ligne_matrix_double(matrix_double *A , unsigned int line1 , unsigned int line2) {
    for(int j=0 ; j<A->nb_col ; j++)
        swap(&A->values[line1*A->nb_col + j] , &A->values[line2*A->nb_col + j]);
}

/* ********************************************************************************************************************** */

void swap_col_matrix_double(matrix_double *A , unsigned int col1 , unsigned int col2) {
    for(int i=0 ; i<A->nb_line ; i++)
        swap(&A->values[i*A->nb_col + col1] , &A->values[i*A->nb_col + col2]);
}