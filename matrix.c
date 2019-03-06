//
// Created by root on 07/12/18.
//

#include "matrix.h"

void init_matrix(matrix *matrix , unsigned int nb_line_p , unsigned int nb_col_p) {
    if(nb_line_p <= 0 || nb_col_p <= 0) {
        perror("Try to create a matrix without line or column");
        return;
    }

    matrix->values = (long*)calloc(sizeof(double) , nb_line_p*nb_col_p);
    for(int i=0 ; i<nb_line_p ; i++) {
        for(int j=0 ; j<nb_col_p ; j++) {
            matrix->values[i*nb_col_p + j] = 0;
        }
    }

    matrix->nb_line = nb_line_p;
    matrix->nb_col = nb_col_p;
}

/* ********************************************************************************************************************** */


void setCoeff_matrix(matrix *matrix , unsigned int line , unsigned int col , long newvalue) {
    if(line < 0 || line >= matrix->nb_line || col < 0 || col >= matrix->nb_col) {
        perror("Matrix index not correct");
        return;
    }
    MATRIX_P(matrix,line,col) = newvalue;
}

/* ********************************************************************************************************************** */

void setAllCoeff_matrix(matrix *matrix , long value) {
    for(unsigned int i=0 ; i<matrix->nb_line ; i++) {
        for(unsigned int j=0 ; j<matrix->nb_col ; j++) {
            //matrix->values[i*matrix->nb_col + j] = value;
            MATRIX_P(matrix,i,j) = value;
        }
    }
}

/* ********************************************************************************************************************** */

/**
 * Set the coefficients of a matrix line by line from an array
 * @param matrix
 * @param coeffs
 * @param sizeArray
 */
void setCoeff_matrix_array(matrix *matrix , long *coeffs, unsigned int sizeArray) {
    if(matrix->nb_line*matrix->nb_col != sizeArray) {
        perror("SetCoeffArray : Error array size not correct");
        exit(EXIT_FAILURE);
    }

    for(unsigned int i=0 ; i<sizeArray ; i++)
        matrix->values[i] = coeffs[i];

}

/* ********************************************************************************************************************** */

void destroy_matrix(matrix matrix) {
    free(matrix.values);
}

/* ********************************************************************************************************************** */

void print_matrix(matrix matrix , char *name) {
    printf("%s :\n",name);
    for(int i=0 ; i<matrix.nb_line ; i++) {
        for(int j=0 ; j<matrix.nb_col ; j++) {
            if(j == 0) printf("| ");

            if(matrix.values[i*matrix.nb_col + j] == 0) printf(" 0.00 ");
            if(matrix.values[i*matrix.nb_col + j] > 0) printf("+");
            if(matrix.values[i*matrix.nb_col + j] != 0) printf("%.2ld ",matrix.values[i*matrix.nb_col + j]);
            if(matrix.values[i*matrix.nb_col + j] < 10  && matrix.values[i*matrix.nb_col + j] > -10) printf(" ");
            if(j == matrix.nb_col - 1) printf("|\n");
        }
    }

    printf("\n");
}

/* ********************************************************************************************************************** */

void change_nb_line_matrix(matrix *matrix , unsigned int new_nb_line) {
    if(new_nb_line == matrix->nb_line) return;

    long *new_values = (long*)calloc(sizeof(double) , new_nb_line*matrix->nb_col);

    if(new_nb_line < matrix->nb_line) {
        for(int i=0 ; i<matrix->nb_line ; i++) {
            for(int j=0 ; j<matrix->nb_col ; j++) {
                if(i < new_nb_line) {
                    //new_values[i*matrix->nb_col + j] = matrix->values[i*matrix->nb_col + j];
                    new_values[i*matrix->nb_col + j] = MATRIX_P(matrix,i,j);
                }
            }
        }
    }
    else {
        for(int i=0 ; i<new_nb_line ; i++) {
            for(int j=0 ; j<matrix->nb_col ; j++) {
                if(i <matrix->nb_line) {
                    //new_values[i*matrix->nb_col + j] = matrix->values[i*matrix->nb_col + j];
                    new_values[i*matrix->nb_col + j] = MATRIX_P(matrix,i,j);
                }
                else {
                    //new_values[i*matrix->nb_col + j] = 0;
                    MATRIX_P(matrix,i,j) = 0;
                }
            }
        }
    }
    free(matrix->values);
    matrix->values = new_values;
    matrix->nb_line = new_nb_line;
}

/* ********************************************************************************************************************** */

void change_nb_col_matrix(matrix *matrix , unsigned int new_nb_col) {
    if(new_nb_col == matrix->nb_col) return;

    long *newvalues = (long*)calloc(sizeof(double) , matrix->nb_line*new_nb_col);

    if(new_nb_col < matrix->nb_col) {
        for(int i=0 ; i<matrix->nb_line ; i++) {
            for(int j=0 ; j<matrix->nb_col ; j++) {
                if(j < new_nb_col) {
                    //newvalues[i*new_nb_col + j] = matrix->values[i*matrix->nb_col + j];
                    newvalues[i*new_nb_col + j] = MATRIX_P(matrix,i,j);
                }
            }
        }
    }
    else {
        for(int i=0 ; i<matrix->nb_line ; i++) {
            for(int j=0 ; j<new_nb_col ; j++) {
                if(j < matrix->nb_col) {
                    //newvalues[i*new_nb_col + j] = matrix->values[i*matrix->nb_col + j];
                    newvalues[i*new_nb_col + j] = MATRIX_P(matrix,i,j);
                }
                else {
                    newvalues[i*new_nb_col + j] = 0;
                }
            }
        }
    }

    free(matrix->values);
    matrix->values = newvalues;
    matrix->nb_col = new_nb_col;
}

/* ********************************************************************************************************************** */

void change_dim_matrix(matrix *matrix , unsigned int new_nb_line , unsigned int new_nb_col) {
    if(new_nb_line != matrix->nb_line) change_nb_line_matrix(matrix , new_nb_line);
    if(new_nb_col != matrix->nb_col)   change_nb_col_matrix(matrix , new_nb_col);
}

/* ********************************************************************************************************************** */

void add_matrix(matrix *res , matrix A , matrix B) {
    if(A.nb_line != B.nb_line || A.nb_col != B.nb_col) {
        perror("Add matrix : error dimensions");
        return;
    }

    change_dim_matrix(res , A.nb_line , A.nb_col);

    for(unsigned int i=0 ; i<res->nb_line ; i++) {
        for(unsigned int j=0 ; j<res->nb_col ; j++)
            MATRIX_P(res,i,j) = MATRIX(A,i,j) + MATRIX(B,i,j);
    }
}

/* ********************************************************************************************************************** */

void sub_matrix(matrix *res , matrix A , matrix B) {
    if(A.nb_line != B.nb_line || A.nb_col != B.nb_col) {
        perror("Sub matrix : error dimensions");
        return;
    }

    change_dim_matrix(res , A.nb_line , A.nb_col);

    matrix B_neg;
    init_matrix(&B_neg , B.nb_line , B.nb_col);

    scalar_mul_matrix(&B_neg , B , -1);
    add_matrix(res , A , B_neg);
    destroy_matrix(B_neg);
}

/* ********************************************************************************************************************** */

void scalar_mul_matrix(matrix *res , matrix A , long lambda) {

    change_dim_matrix(res, A.nb_line, A.nb_col);

    for (unsigned int i = 0; i < res->nb_line; i++) {
        for (unsigned int j = 0; j < res->nb_col; j++)
            MATRIX_P(res, i, j) = lambda * MATRIX(A, i, j);
    }
}

/* ********************************************************************************************************************** */

void scalar_div_matrix(matrix *res , matrix A , long lambda) {
    if(!lambda) {
        perror("Scalar div : lambda must be != 0");
        return;
    }
    scalar_mul_matrix(res,A,1/lambda);
}

/* ********************************************************************************************************************** */

void dot_product(long *dot , matrix u , matrix v) {
    *dot = 0;

    if(u.nb_line != 1 || v.nb_col != 1 || u.nb_col != v.nb_line) {
        perror("Error dot product : dimensions not respected");
        return;
    }

    for(unsigned int i=0 ; i<u.nb_col ; i++)
        *dot += u.values[i]*v.values[i];
}

/* ********************************************************************************************************************** */

void mul_matrix(matrix *res , matrix A , matrix B) {
    if(A.nb_col != B.nb_line) {
        perror("Error dimension mul_matrix");
        exit(EXIT_FAILURE);
    }

    change_dim_matrix(res , A.nb_line , B.nb_col);

    for(int i=0 ; i<res->nb_line ; i++) {
        for(int j=0 ; j<res->nb_col ; j++) {

            //res->values[i*res->nb_col + j] = 0;
            MATRIX_P(res,i,j) = 0;

            for(int k=0 ; k<A.nb_line ; k++) {
                //res->values[i*res->nb_col+j] += A.values[i*A.nb_col + k]*B.values[k*B.nb_col + j];
                MATRIX_P(res,i,j) += MATRIX(A,i,k)*MATRIX(B,k,j);
            }
        }
    }
}

/* ********************************************************************************************************************** */

void copy_matrix(matrix *DEST , matrix SRC) {
    change_dim_matrix(DEST , SRC.nb_line , SRC.nb_col);

    for(int i=0 ; i<SRC.nb_line ; i++) {
        for (int j = 0; j < SRC.nb_col; j++) {
            //DEST->values[i*SRC.nb_col + j] = SRC.values[i*SRC.nb_col + j];
            MATRIX_P(DEST,i,j) = MATRIX(SRC,i,j);
        }
    }
}

/* ********************************************************************************************************************** */

void identity_matrix(matrix *id, unsigned int size) {
    change_dim_matrix(id, size, size);

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

void matrix_line_permutation(matrix *P , unsigned int size , unsigned int line1 , unsigned int line2) {
    identity_matrix(P , size);
    swap_ligne_matrix(P , line1 , line2);
}

/* ********************************************************************************************************************** */

void matrix_col_permutation(matrix *Q ,unsigned int size ,  unsigned int col1 , unsigned int col2) {
    identity_matrix(Q , size);
    swap_col_matrix(Q , col1 , col2);
}

/* ********************************************************************************************************************** */

void transpose_matrix(matrix *transpose , matrix A) {
    change_dim_matrix(transpose , A.nb_col , A.nb_line);

    for(unsigned int i=0 ; i< A.nb_line ; i++) {
        for(unsigned int j=0 ; j<A.nb_col ; j++) {
            setCoeff_matrix(transpose , j , i , A.values[i*A.nb_col + j]);
        }
    }
}


/* ********************************************************************************************************************** */

void LU_decomposition_matrix(matrix *L , matrix *U , matrix A) {

    if(L != NULL)
        identity_matrix(L , A.nb_line); // L = Id
    copy_matrix(U , A);             // U = A

    long coeff, temp;

    for(int j=0 ; j<A.nb_col ; j++) {
        for(int i=j+1 ; i<A.nb_line ; i++) {

            coeff = -1 * MATRIX_P(U,i,j) / MATRIX_P(U,j,j);

            for(int k=j ; k<A.nb_col ; k++) {
                temp = MATRIX_P(U,j,k) * coeff;

                if(L != NULL) {
                    MATRIX_P(L,i,j) = -coeff;
                }

                MATRIX_P(U,i,k) += temp;
            }
        }
    }
}

/* ********************************************************************************************************************** */

/**
 * Swap two lines
 * Index start from 0
 * @param A
 * @param line1
 * @param line2
 */
void swap_ligne_matrix(matrix *A , unsigned int line1 , unsigned int line2) {
    if(line1 != line2) {
        for(int j=0 ; j<A->nb_col ; j++)
            swap_long(&A->values[line1*A->nb_col + j] , &A->values[line2*A->nb_col + j]);
    }

}

/* ********************************************************************************************************************** */

/**
 * Swap two columns
 * Index start from 0
 * @param A
 * @param col1
 * @param col2
 */
void swap_col_matrix(matrix *A , unsigned int col1 , unsigned int col2) {
    if(col1 != col2) {
        for(int i=0 ; i<A->nb_line ; i++)
            swap_long(&A->values[i*A->nb_col + col1] , &A->values[i*A->nb_col + col2]);
    }

}

/* ********************************************************************************************************************** */

/**
 * Return the max element in a submatrix of A
 * Set the index of the max in max_index_line and max_index_col
 * @param max_index_line
 * @param max_index_col
 * @param A
 * @param submatrix_index_line
 * @param submatrix_index_col
 * @return
 */
long index_max_submatrix(unsigned int *max_index_line , unsigned int *max_index_col , matrix A , unsigned int submatrix_index_line , unsigned int submatrix_index_col) {
    double max = MATRIX(A,submatrix_index_line,submatrix_index_col);
    *max_index_line = submatrix_index_line;
    *max_index_col = submatrix_index_col;

    for (unsigned int i=submatrix_index_line ; i<A.nb_line ; i++) {
        for(unsigned int j=submatrix_index_col ; j<A.nb_col ; j++) {
            //if(A.values[i*A.nb_col + j] > max) {
            if(MATRIX(A,i,j) > max) {
                max = MATRIX(A,i,j);
                *max_index_line = i;
                *max_index_col = j;
            }
        }
    }
    return max;
}

/* ********************************************************************************************************************** */

/**
 * Input : a matrix of size m*1
 * Output : the norm -> sqrt(sum(xi^2))
 * @param matrix
 * @return
 */
double norm_vector(matrix matrix) {
    if(matrix.nb_col != 1) {
        perror("Norm must be called on vectors only");
        return -1;
    }

    double norm = 0;

    for(unsigned int i=0 ; i<matrix.nb_line ; i++)
        norm += pow(matrix.values[i] , 2);

    norm = sqrt(norm);
    return norm;
}

/* ********************************************************************************************************************** */

/**
 * Replace one column of a matrix with another vector
 * @param A
 * @param vector
 * @param indexColumn
 */
void setColumn_matrix(matrix *A , matrix vector , unsigned int indexColumn) {
    if(vector.nb_col != 1 || vector.nb_line != A->nb_line) {
        perror("Set column : second argument must be a vector with correct dimensions");
        return;
    }

    if(indexColumn < 0 || indexColumn >= A->nb_col) {
        perror("Set column : index is not valid");
        return;
    }

    for(unsigned int i=0 ; i<vector.nb_line ; i++)
        MATRIX_P(A,i,indexColumn) = vector.values[i];

}

/* ********************************************************************************************************************** */

/**
 * Extract one column from a matrix
 * @param column
 * @param A
 * @param indexColumn
 */
void getColum_matrix(matrix *column , matrix A , unsigned int indexColumn) {
    if(indexColumn < 0 || indexColumn >= A.nb_col) {
        perror("getColumn : index is not valid");
        return;
    }

    change_dim_matrix(column , A.nb_line , 1);

    for(unsigned int i=0 ; i<A.nb_line ; i++)
        column->values[i] = MATRIX(A,i,indexColumn);
}

/* ********************************************************************************************************************** */

