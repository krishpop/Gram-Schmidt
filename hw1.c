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
void double_gs(double *a, int n, double *u, double *g);
void transpose(double *matrix, double **transposed_matrix, int n);
void multiply(double *m1, double *m2, double **new_matrix, int n);
double *projection(double dp, double *m1, double *m2, int v1, int v2, int n);
double dot_product(double *m1, double *m2, int v1, int v2, int n);
void normalize(double v_mag, double *m, int v, int n);
double vector_magn(double *m, int v, int n);
void print_matrix(double *m, int rows, int columns);
void identity_matrix(double** m, int n);

void inv_double_gs(double *a, int n, double *u, double *b) {
    double *g;
    identity_matrix(&g, n);
    double_gs(a, n, u, g);
    double *transp_u;
    transpose(u, &transp_u, n);
    multiply(g, transp_u, &b, n);
}

void double_gs(double *a, int n, double *u, double *g) {
    int v_i, u_j, row, double_gs;
    double *v_matrix, *proj_uv, *proj_ui;
    double dp;
    for (v_i=0; v_i<n; v_i++) {
        for (row=0; row<n; row++) {
            u[row*n+v_i]=a[row*n+v_i];
        }
        double_gs = 0;
        v_matrix = a;
        for (u_j=v_i-1; u_j>=0; u_j--) {
            dp = dot_product(v_matrix, u, v_i, u_j, n);
            proj_uv = projection(dp, v_matrix, u, v_i, u_j, n);
            for (row=0; row<n; row++) {
                u[row*n+v_i] = u[row*n+v_i] - proj_uv[row];
                g[row*n+v_i] = g[row*n+v_i] - dp*g[row*n+u_j];
            }
            free(proj_uv);
            if (u_j == 0 && double_gs == 0) {
                u_j = v_i-1;
                v_matrix = u;
                double_gs = 1;
            }
        }
        // normalize the vector in g according to magnitude of u_v
        double v_mag = vector_magn(u, v_i, n);
        for (row=0; row<n; row++) {
            g[row*n+v_i] = g[row*n+v_i] / v_mag;
        }
        normalize(v_mag, u, v_i, n);
    }
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

double *projection(double dpsum, double *a, double *u, int v_i, int u_j, int n){
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

void normalize(double v_mag, double *u, int v, int n) {
    int row;
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

int main() {
    double *a, *b, *u;
    int n = 3;
    a = (double*)malloc(sizeof(double)*n*n);
    b = (double*)malloc(sizeof(double)*n*n);
    u = (double*)malloc(sizeof(double)*n*n);
    a[0] = 1;
    a[1] = 1;
    a[2] = 1;
    a[3] = -1;
    a[4] = 0;
    a[5] = 1;
    a[6] = 1;
    a[7] = 1;
    a[8] = 2;
    print_matrix(a, n, n);
    inv_double_gs(a, n, u, b);
    print_matrix(b, n, n);
    return 0;
}
