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
    //matrixDouble->values[line*matrixDouble->nb_col + col] = newvalue;
    MATRIX_P(matrixDouble,line,col) = newvalue;
}

/* ********************************************************************************************************************** */

void setAllCoeff_matrix_double(matrix_double *matrixDouble , double value) {
    for(unsigned int i=0 ; i<matrixDouble->nb_line ; i++) {
        for(unsigned int j=0 ; j<matrixDouble->nb_col ; j++) {
            //matrixDouble->values[i*matrixDouble->nb_col + j] = value;
            MATRIX_P(matrixDouble,i,j) = value;
        }
    }
}

/* ********************************************************************************************************************** */

/**
 * Set the coefficients of a matrix line by line from an array
 * @param matrixDouble
 * @param coeffs
 * @param sizeArray
 */
void setCoeff_matrix_double_array(matrix_double *matrixDouble , double *coeffs, unsigned int sizeArray) {
    if(matrixDouble->nb_line*matrixDouble->nb_col != sizeArray) {
        perror("SetCoeffArray : Error array size not correct");
        exit(EXIT_FAILURE);
    }

    for(unsigned int i=0 ; i<sizeArray ; i++)
        matrixDouble->values[i] = coeffs[i];

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
                    //new_values[i*matrixDouble->nb_col + j] = matrixDouble->values[i*matrixDouble->nb_col + j];
                    new_values[i*matrixDouble->nb_col + j] = MATRIX_P(matrixDouble,i,j);
                }
            }
        }
    }
    else {
        for(int i=0 ; i<new_nb_line ; i++) {
            for(int j=0 ; j<matrixDouble->nb_col ; j++) {
                if(i <matrixDouble->nb_line) {
                    //new_values[i*matrixDouble->nb_col + j] = matrixDouble->values[i*matrixDouble->nb_col + j];
                    new_values[i*matrixDouble->nb_col + j] = MATRIX_P(matrixDouble,i,j);
                }
                else {
                    //new_values[i*matrixDouble->nb_col + j] = 0;
                    MATRIX_P(matrixDouble,i,j) = 0;
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
                    //newvalues[i*new_nb_col + j] = matrixDouble->values[i*matrixDouble->nb_col + j];
                    newvalues[i*new_nb_col + j] = MATRIX_P(matrixDouble,i,j);
                }
            }
        }
    }
    else {
        for(int i=0 ; i<matrixDouble->nb_line ; i++) {
            for(int j=0 ; j<new_nb_col ; j++) {
                if(j < matrixDouble->nb_col) {
                    //newvalues[i*new_nb_col + j] = matrixDouble->values[i*matrixDouble->nb_col + j];
                    newvalues[i*new_nb_col + j] = MATRIX_P(matrixDouble,i,j);
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

void add_matrix_double(matrix_double *res , matrix_double A , matrix_double B) {
    if(A.nb_line != B.nb_line || A.nb_col != B.nb_col) {
        perror("Add matrix : error dimensions");
        return;
    }

    change_dim_matrix_double(res , A.nb_line , A.nb_col);

    for(unsigned int i=0 ; i<res->nb_line ; i++) {
        for(unsigned int j=0 ; j<res->nb_col ; j++)
            MATRIX_P(res,i,j) = MATRIX(A,i,j) + MATRIX(B,i,j);
    }
}

/* ********************************************************************************************************************** */

void sub_matrix_double(matrix_double *res , matrix_double A , matrix_double B) {
    if(A.nb_line != B.nb_line || A.nb_col != B.nb_col) {
        perror("Sub matrix : error dimensions");
        return;
    }

    change_dim_matrix_double(res , A.nb_line , A.nb_col);

    matrix_double B_neg;
    init_matrix_double(&B_neg , B.nb_line , B.nb_col);

    scalar_mul_matrix_double(&B_neg , B , -1);
    add_matrix_double(res , A , B_neg);
    destroy_matrix_double(B_neg);
}

/* ********************************************************************************************************************** */

void scalar_mul_matrix_double(matrix_double *res , matrix_double A , double lambda) {

    change_dim_matrix_double(res, A.nb_line, A.nb_col);

    for (unsigned int i = 0; i < res->nb_line; i++) {
        for (unsigned int j = 0; j < res->nb_col; j++)
            MATRIX_P(res, i, j) = lambda * MATRIX(A, i, j);
    }
}

/* ********************************************************************************************************************** */

void scalar_div_matrix_double(matrix_double *res , matrix_double A , double lambda) {
    if(!lambda) {
        perror("Scalar div : lambda must be != 0");
        return;
    }
    scalar_mul_matrix_double(res,A,1/lambda);
}

/* ********************************************************************************************************************** */

void dot_product(double *dot , matrix_double u , matrix_double v) {
    *dot = 0;

    if(u.nb_line != 1 || v.nb_col != 1 || u.nb_col != v.nb_line) {
        perror("Error dot product : dimensions not respected");
        return;
    }

    for(unsigned int i=0 ; i<u.nb_col ; i++)
        *dot += u.values[i]*v.values[i];
}

/* ********************************************************************************************************************** */

void mul_matrix_double(matrix_double *res , matrix_double A , matrix_double B) {
    if(A.nb_col != B.nb_line) {
        perror("Error dimension mul_matrix_double");
        exit(EXIT_FAILURE);
    }

    change_dim_matrix_double(res , A.nb_line , B.nb_col);

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

void copy_matrix_double(matrix_double *DEST , matrix_double SRC) {
    change_dim_matrix_double(DEST , SRC.nb_line , SRC.nb_col);

    for(int i=0 ; i<SRC.nb_line ; i++) {
        for (int j = 0; j < SRC.nb_col; j++) {
            //DEST->values[i*SRC.nb_col + j] = SRC.values[i*SRC.nb_col + j];
            MATRIX_P(DEST,i,j) = MATRIX(SRC,i,j);
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

void matrix_line_permutation(matrix_double *P , unsigned int size , unsigned int line1 , unsigned int line2) {
    identity_matrix_double(P , size);
    swap_ligne_matrix_double(P , line1 , line2);
}

/* ********************************************************************************************************************** */

void matrix_col_permutation(matrix_double *Q ,unsigned int size ,  unsigned int col1 , unsigned int col2) {
    identity_matrix_double(Q , size);
    swap_col_matrix_double(Q , col1 , col2);
}

/* ********************************************************************************************************************** */

void transpose_matrix_double(matrix_double *transpose , matrix_double A) {
    change_dim_matrix_double(transpose , A.nb_line , A.nb_col);

    for(unsigned int i=0 ; i< A.nb_line ; i++) {
        for(unsigned int j=0 ; j<A.nb_col ; j++) {
            setCoeff_matrix_double(transpose , j , i , A.values[i*A.nb_col + j]);
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
void swap_ligne_matrix_double(matrix_double *A , unsigned int line1 , unsigned int line2) {
    if(line1 != line2) {
        for(int j=0 ; j<A->nb_col ; j++)
            swap(&A->values[line1*A->nb_col + j] , &A->values[line2*A->nb_col + j]);
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
void swap_col_matrix_double(matrix_double *A , unsigned int col1 , unsigned int col2) {
    if(col1 != col2) {
        for(int i=0 ; i<A->nb_line ; i++)
            swap(&A->values[i*A->nb_col + col1] , &A->values[i*A->nb_col + col2]);
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
double index_max_submatrix_double(unsigned int *max_index_line , unsigned int *max_index_col , matrix_double A , unsigned int submatrix_index_line , unsigned int submatrix_index_col) {
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
 * @param matrixDouble
 * @return
 */
double norm_vector_double(matrix_double matrixDouble) {
    if(matrixDouble.nb_col != 1) {
        perror("Norm must be called on vectors only");
        return -1;
    }

    double norm = 0;

    for(unsigned int i=0 ; i<matrixDouble.nb_line ; i++)
        norm += pow(matrixDouble.values[i] , 2);

    norm = sqrt(norm);
    return norm;
}

/* ********************************************************************************************************************** */

void PLUQ_decomposition(matrix_double *P , matrix_double *L , matrix_double *U , matrix_double *Q , matrix_double A) {
    double coeff;
    unsigned int size = A.nb_line;
    unsigned int max_index_line, max_index_col;

    matrix_double R,S,temp;
    init_matrix_double(&R , A.nb_line , A.nb_col);
    init_matrix_double(&S , A.nb_line , A.nb_col);
    init_matrix_double(&temp , A.nb_line , A.nb_col);

    identity_matrix_double(L , size);
    identity_matrix_double(P , size);
    identity_matrix_double(Q , size);
    copy_matrix_double(U , A);

    for(unsigned int j=0 ; j<U->nb_col ; j++) {
        index_max_submatrix_double(&max_index_line , &max_index_col , *U , j , j);

        swap_ligne_matrix_double(U , j , max_index_line);
        matrix_line_permutation(&R , size , j , max_index_line);

        swap_col_matrix_double(U , j , max_index_col);
        matrix_col_permutation(&S , size , j , max_index_col);

        mul_matrix_double(&temp , S , *Q);
        copy_matrix_double(Q , temp);

        mul_matrix_double(&temp , *L , R);
        copy_matrix_double(L , temp);

        mul_matrix_double(&temp, R , *L);
        copy_matrix_double(L , temp);

        mul_matrix_double(&temp , *P , R);
        copy_matrix_double(P , temp);


        for(int i=j+1 ; i<U->nb_line ; i++) {
            coeff = MATRIX_P(U,i,j) / MATRIX_P(U,j,j);
            for (int k = j; k < U->nb_col; k++) MATRIX_P(U,i,k) -= coeff * MATRIX_P(U,j,k);
            MATRIX_P(L,i,j) = coeff;
        }
    }

    destroy_matrix_double(R);
    destroy_matrix_double(S);
    destroy_matrix_double(temp);
}

/* ********************************************************************************************************************** */

/**
 * Replace one column of a matrix with another vector
 * @param A
 * @param vector
 * @param indexColumn
 */
void setColumn_matrix_double(matrix_double *A , matrix_double vector , unsigned int indexColumn) {
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
void getColum_matrix_double(matrix_double *column , matrix_double A , unsigned int indexColumn) {
    if(indexColumn < 0 || indexColumn >= A.nb_col) {
        perror("getColumn : index is not valid");
        return;
    }

    change_dim_matrix_double(column , A.nb_line , 1);

    for(unsigned int i=0 ; i<A.nb_line ; i++)
        column->values[i] = MATRIX(A,i,indexColumn);
}

/* ********************************************************************************************************************** */

/**
 * Input : A of size m*n
 * Output : QR factorization
 * @param Q
 * @param R
 * @param A
 */
void QR_Givens(matrix_double *Q , matrix_double *R , matrix_double A) {
    unsigned int m = A.nb_line , n = A.nb_col;
    double c,s;

    matrix_double G,Gt,temp;
    init_matrix_double(&temp , m , n);
    init_matrix_double(&G , m , m);
    init_matrix_double(&Gt , m ,m);

    identity_matrix_double(Q , A.nb_line);
    copy_matrix_double(R , A);

    for(unsigned int j=0 ; j<n ; j++) {
        for(unsigned int i=j+1 ; i<m ; i++) {

            /* ******************************************************************** */
            // Build G

            setAllCoeff_matrix_double(&G , 0);

            for(unsigned int k=0 ; k<m ; k++) {
                if(k != i && k != j) setCoeff_matrix_double(&G , k , k , 1);
            }

            c = MATRIX_P(R,j,j) / (sqrt(pow(MATRIX_P(R,j,j) , 2) + pow(MATRIX_P(R,i,j) , 2)));
            s = MATRIX_P(R,i,j) / (sqrt(pow(MATRIX_P(R,j,j) , 2) + pow(MATRIX_P(R,i,j) , 2)));

            setCoeff_matrix_double(&G , i , i , c);
            setCoeff_matrix_double(&G , j , j , c);

            setCoeff_matrix_double(&G , i , j , -s);
            setCoeff_matrix_double(&G , j , i , s);

            /* ******************************************************************** */

            // update R
            mul_matrix_double(&temp , G , *R);
            copy_matrix_double(R , temp);

            // update Q
            transpose_matrix_double(&Gt , G);
            mul_matrix_double(&temp , *Q , Gt);
            copy_matrix_double(Q , temp);
        }
    }

    destroy_matrix_double(G);
    destroy_matrix_double(Gt);
    destroy_matrix_double(temp);
}

/* ********************************************************************************************************************** */

/**
 * Input : A of size m*n
 * Output : QR factorization
 * @param Q
 * @param R
 * @param A
 */
void QR_Gram_Schimdt(matrix_double *Q , matrix_double *R , matrix_double A) {

    unsigned int m = A.nb_line , n= A.nb_col;

    change_dim_matrix_double(Q , A.nb_line , A.nb_line);

    change_dim_matrix_double(R , m , n);
    setAllCoeff_matrix_double(R , 0);

    matrix_double a,q_j,q_i,q_iT;
    init_matrix_double(&a , A.nb_line , 1);
    init_matrix_double(&q_j , A.nb_line , 1);
    init_matrix_double(&q_i , A.nb_line , 1);
    init_matrix_double(&q_iT , 1 , A.nb_line);


    getColum_matrix_double(&a , A , 0);
    MATRIX_P(R,0,0) = norm_vector_double(a);

    scalar_div_matrix_double(&q_j , a , MATRIX_P(R,0,0));
    setColumn_matrix_double(Q , q_j , 0);


    for(unsigned int j=1 ; j<A.nb_col ; j++) {

        // q(j) = a(j)
        getColum_matrix_double(&q_j , A , j);

        for(unsigned int i=0 ; i<j ; i++) {
            getColum_matrix_double(&q_i, *Q, i);
            transpose_matrix_double(&q_iT, q_i);

            dot_product(&MATRIX_P(R, i, j), q_iT, q_j);

            /* **************************************** */
            // TODO
            // q(j) <- q(j) - R(i,j)*q(j)
            /* **************************************** */


        }

        MATRIX_P(R,j,j) = norm_vector_double(q_j);
        scalar_div_matrix_double(&q_j , MATRIX_P(R,j,j));

        setColumn_matrix_double(Q, q_j, j);
    }

    
    destroy_matrix_double(a);
    destroy_matrix_double(q_j);
    destroy_matrix_double(q_i);
    destroy_matrix_double(q_iT);

}