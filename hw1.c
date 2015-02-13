//
//  CPSC 440 - hw1.c
//  Double Gram-Schmidt Inversion
//
//  Created by Krishnan Srinivasan on 2/8/15.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void inv_double_gs(double *a, int n, double *u, double *b);
void transpose(double *matrix, double **transposed_matrix, int n);
void multiply(double *m1, double *m2, double **new_matrix, int n);
double *projection(double *m1, double *m2, int v1, int v2, int n);
double dot_product(double *m1, double *m2, int v1, int v2, int n);
void normalize(double *m, int v, int n);
double vector_magn(double *m, int v, int n);
void print_matrix(double *m, int rows, int columns);
void identity_matrix(double** m, int n);

void inv_double_gs(double *a, int n, double *u, double *b) {
    u = (double*)malloc(sizeof(double)*n*n);
    int v_i, u_j;
    for (v_i=0; v_i<n; v_i++) {
        int row;
        for (row=0; row<n; row++) {
            u[row*n+v_i]=a[row*n+v_i];
        }
        int double_gs = 0;
        double *a_gs = a;
        for (u_j=v_i-1; u_j>=0; u_j--) {
            double *proj_uv = projection(a_gs, u, v_i, u_j, n);
            for (row=0; row<n; row++) {
                u[row*n+v_i]-=proj_uv[row]; 
            }
            free(proj_uv);
            if (u_j == 0 && double_gs == 0) {
                u_j = v_i-1;
                a_gs = u;
                double_gs = 1;
            }
        }
        normalize(u, v_i, n);
    }
    print_matrix(u, n, n);
}

void transpose(double *a, double **transposed_matrix, int n) {
    int row, column;
    *transposed_matrix = (double*)malloc(sizeof(double)*n*n);
    
    for (row=0; row < n; row++) {
        for (column=0; column < n; column++) {
            (*transposed_matrix)[row*n+column] = a[column*n+row]; 
            (*transposed_matrix)[column*n+row] = a[row*n+column];
        }
    }
}

double *projection(double *a, double *u, int v_i, int u_j, int n) {
    double dpsum = dot_product(a, u, v_i, u_j, n);
    int row;
    double *proj = malloc(sizeof(double)*n);
    for (row=0; row<n; row++) {
        proj[row] = u[row*n+u_j]*dpsum;
    }
    return proj;
}

double dot_product(double *a, double *b, int v1, int v2, int n) {
    int row;
    double sum = 0;
    for (row=0; row<n; row++) {
        sum += a[row*n + v1] * b[row*n + v2];
    }
    return sum;
}

void normalize(double *u, int v, int n) {
    int row;
    double v_mag = vector_magn(u, v, n);
    for (row=0; row<n; row++) {
        u[row*n+v] = u[row*n+v]/v_mag;
    }
}

double vector_magn(double *a, int v, int n) {
    int row;
    double sum_squared = dot_product(a, a, v, v, n);
    return sqrt(sum_squared);
}

void multiply(double *m1, double *m2, double **new_matrix, int n) {
    int new_row, new_column, old_x;
    *new_matrix = (double*)malloc(sizeof(double)*n*n);

    for (new_row = 0; new_row < n; new_row++) {
        
        for (new_column = 0; new_column < n; new_column++) {
            (*new_matrix)[new_row*n+new_column] = 0;
            
            for (old_x = 0; old_x < n; old_x++) {
                (*new_matrix)[new_row*n+new_column] += m1[new_row*n+old_x] \
                    * m2[old_x*n+new_column];
            }
        }
    }
}

void identity_matrix(double **m, int n) {
    int i;
    *m = (double*)calloc(sizeof(double), n*n);
    for (i=0; i<n; i++) {
        (*m)[i*n+i] = (double)1;
    }
}

int main() {
    int n = 3;
    double *a = malloc(sizeof(double)*n*n);
    double *u, *b;
    int i, j;
    a[0] = a[6] = a[1] = a[7] = a[2] = a[5] = 1;
    a[3] = -1;
    a[4] = 0;
    a[8] = 2;
    print_matrix(a, n, n);
    inv_double_gs(a, n, u, b);    
    return 0;
}

void print_matrix(double *m, int rows, int columns){
    int row, column;
    for (row=0; row<rows; row++) {
        for (column=0; column<columns; column++) {
            printf("[%d, %d]: %f\t", row, column, m[row*rows+column]);
        }
        printf("\n");
    }
    printf("\n");
}

