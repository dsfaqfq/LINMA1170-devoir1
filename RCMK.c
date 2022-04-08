#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {
    double **a;
    double *data;
    int n, k1, k2;
    bool LU;
} band_matrix;

typedef struct node {
    struct node * next;
    double value;
} node_t;

typedef struct queue {
    struct node* tail;
    int size;
} queue_t;

typedef struct {
    int m, n; // dimensions de la matrice
    double* data; // tableau 1D de taille m*n contenant les entrées de la matrice
    double** a; // tableau 1D de m pointeurs vers chaque ligne, pour pouvoir appeler a[i][j]
    bool LU;
} matrix;

matrix* allocate_matrix(int m, int n) {
    matrix* mat = (matrix*)malloc(sizeof(matrix));
    mat->m = m, mat->n = n;
    mat->data = (double*)malloc(m * n * sizeof(double));
    if (mat->data == NULL) return NULL;
    mat->a = (double**)malloc(m * sizeof(double*));
    if (mat->a == NULL) return NULL;
    for (int i = 0; i < m; i++)
        mat->a[i] = mat->data + i * n;
    mat->LU = false;
    return mat;
}

void free_matrix(matrix* mat) {
    if (mat == NULL) return;
    free(mat->a);
    free(mat->data);
}

int* getdegrees(matrix* mat) {
    int n = mat->n;
    double **a = mat->a;
    int* degrees = calloc(mat->n, sizeof(int));
    for (int i = 0 ; i < n ; i++) {
        int degree = 0;
        for (int j = 0 ; j<n ; j++) {
            if (a[i][j] != 0.0) degree ++;
        }
        degrees[i] = degree;
    }
    return degrees;
}

int enqueue(queue_t* q, double val) {
    struct node* new_node = (struct node*) malloc(sizeof(struct node)) ;
    if (new_node == NULL) {
        return -1 ;
    }
    new_node->value = val ;
    if (q->size == 0) {
        new_node->next = new_node ;
        q->tail = new_node ;
    }
    else {
        new_node->next = q->tail->next ;
        q->tail->next = new_node ;
    }
    q->size ++ ;
    return 0 ;
}

double dequeue(queue_t* q) {
    struct node* removed = q->tail ;
    if (q->size == 1) {
        q->tail = NULL ;
    }
    else {
        struct node* current = q->tail ;
        for (int i = 0 ; i<q->size-1; i++){
            current = current->next ;
        }
        current->next = q->tail->next ;
        q->tail = current ;
    }
    double value = removed->value ;
    free(removed) ;
    q->size -- ;
    return value ;
}


band_matrix* allocate_band_matrix (int n, int k1, int k2){
    band_matrix * mat =  malloc(sizeof(band_matrix));
    double *data = calloc (n * (k1 + k2 + 1), sizeof(double));
    double **a = malloc(n * sizeof(double *));
    for (int i = 0;i<n;i++){
        a[i] = data + i * (k1 + k2 + 1);
    }
    mat->a = a;
    mat->data = data;
    mat->k1 = k1;
    mat->k2 = k2;
    mat->n = n;
    mat->LU = false;
    return mat;
}

int LU_factorisation(size_t n, band_matrix* A) {
    if (A->LU) return 0;

    double** a = A->a;
    int k1 = A->k1;
    int k2 = A->k2;

    for (int k = 0; k < n; k++) {
        double akk = a[k][k1];
        if (akk == 0) return -1;
        for (int i = k + 1; i <= k + k1 && i < n; i++) {
            a[i][k1 + k - i] /= akk;
            for (int j = k + 1; j <= k + k2 && j < n; j++) {
                a[i][k1 + j - i] -= a[i][k1 + k - i] * a[k][k1 + j - k];
            }
        }
    }
    A->LU = true;
    return 0;
}

int Back_substitution(int n, band_matrix* U, double* x, double* y) {
    int k1 = U->k1;
    double** u = U->a;
    for (int k = n - 1; k >= 0; k--) {
        x[k] = y[k];
        for (int i = k + 1; i <= k + k1 && i < n; i++) {
            x[k] -= u[k][k1 + i - k] * x[i];
        }
        x[k] /= u[k][k1];
    }
    return 0;
}

int Forward_substitution(int n, band_matrix* L, double* y, double* b) {
    int k1 = L->k1;
    int k2 = L->k2;
    double** l = L->a;
    for (int k = 0; k < n; k++) {
        y[k] = b[k];
        for (int i = k-1 ; i >= k - k1 && i >= 0; i--) {
            y[k] -= l[k][k1 + i - k] * y[i];
        }
    }
    return 0;
}

int solve_system(int n, band_matrix* A, double* b, double* x) {
    double* y = malloc(sizeof(double) * n);
    LU_factorisation(n, A);
    Forward_substitution(n, A, y, b);
    Back_substitution(n, A, x, y);
    free(y);
    return 0;
}

void free_band_matrix (band_matrix *mat){
    free(mat->data);
    free(mat->a);
    free(mat);
}

int main(){
    int k1 = 1;
    int k2 = 1;
    int n = 10;
    band_matrix* A = allocate_band_matrix(n, k1, k2);
    
    for (int i = 0; i < n; i++) {
        A->a[i][k1] = 2.0;
    }
    for (int i = 1; i < n; i++) {
        A->a[i][k1-1] = -1.0;
    }
    for (int i = 0; i < n-1; i++) {
        A->a[i][k1+1] = -1.0;
    }


    double* b = calloc(n, sizeof(double));
    double* x = calloc(n, sizeof(double));

    b[0]   = 1.0;
    b[n-1] = 1.0;

    solve_system(n, A, b, x);

    for (int i = 0 ; i < n; i++) {
        printf("%f ", x[i]);
    }
    printf("\n");

    free_band_matrix(A);
    free(b);
    free(x);

    return 0;
}