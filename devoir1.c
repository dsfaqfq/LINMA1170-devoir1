#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

double L;
double H;
double Lz;
double f; // can change
int M; // can change (from 19 to 99)
int N;
double Mt; // can change (from 8 to 90)
double S11[2];
double S12[2];
double S21[2];
double S22[2];
double* S[] = { S11, S12, S21, S22 };
double dx;
double dy;
double dt;
double mu0;
double sigma;
double R;
int m;

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

int LU_decomposition_pivot(int n, matrix* A, int* p) {
    if (A->LU) return -1;
    for (int k = 0; k < A->m; k++) {

        // choix du pivot
        int idx = k; // indice du pivot
        double max = fabs(A->a[k][k]); //valeur au pivot
        for (int i = k + 1; i < n; i++) {
            if (fabs(A->a[i][k]) > max) { //changement de pivot si valeur observée plus grande
                max = A->a[i][k];
                idx = i;
            }
        }
        if (max == 0) return -1;

        // application du pivot une fois choisi
        double* tmp1 = A->a[k];
        A->a[k] = A->a[idx];
        A->a[idx] = tmp1;

        int tmp2 = p[k];
        p[k] = p[idx];
        p[idx] = tmp2;

        // factorisation LU avec le nouveau pivot
        double akk = A->a[k][k];
        if (akk == 0) return -1;
        for (int i = k + 1; i < n; i++) {
            A->a[i][k] /= akk;
            for (int j = k + 1; j < n; j++) {
                A->a[i][j] -= A->a[i][k] * A->a[k][j];
            }
        }
    }
    A->LU = true;
    return 0;
}

int Forward_substitution_pivot(int n, matrix* A, int* p, double* y, double* b) {
    if (!A->LU) return -1;
    for (int k = 0; k < n; k++) {
        if (A->a[k][k] == 0) return -1;
        y[k] = b[p[k]]; // recuperation du bon indice dans p contenant les pivots effectues dans LU
        for (int i = k - 1; i >= 0; i--) {
            y[k] -= A->a[k][i] * y[i];
        }
    }
    return 0;
}

int Back_substitution_pivot(int n, matrix* A, int* p, double* x, double* y) {
    if (!A->LU) return -1;
    for (int k = n - 1; k >= 0; k--) {
        if (A->a[k][k] == 0) return -1;
        x[k] = y[k];
        for (int i = k + 1; i < n; i++) {
            x[k] -= A->a[k][i] * x[i];
        }
        x[k] /= A->a[k][k];
    }
    return 0;
}

double E1(double t, double f) {
    // tension initiale
    return 0.001 * sin(2*M_PI*f*t);
}

void initialise_A(matrix* A) {
    double** a = A->a;

    /* coefficients lies au laplacien et à j

       1/mu0 * (d2ai,j/dx2 + d2ai,j/dy2) + j = 0
       1/mu0 * (-4*ai,j + ai+1,j + ai-1,j + ai,j+1 + ai,j-1)/dx2 + j = 0

       j = 0 en dehors de spires
    */
    double coeff = 1.0 / (dx*dx*mu0);
    for (int i = 0; i < M * N; i++) { // pour le laplacien
        a[i][i] = -4.0 * coeff;
        // serie de "if" pour verifier si on est sur un cote/coin du maillage
        if (i % M != 0) a[i][i - 1] = coeff;
        if (i % M != M - 1) a[i][i + 1] = coeff; 
        if (i / M != 0) a[i][i - M] = coeff;
        if (i / M != M - 1) a[i][i + M] = coeff;
    }
    
    double coeff1 = -sigma / dt;
    double coeff2 = sigma  / Lz;
    for (int i = 0; i < 4; i++) { // pour j
        int minx = ((int)(S[i][0] / dx)) - 1;
        int miny = ((int)(S[i][1] / dy)) - 1;
        int maxx = minx + ((int)(0.02 / dx));
        int maxy = miny + ((int)(0.01 / dx));
        for (int j = minx; j <= maxx; j++) {
            for (int k = miny; k <= maxy; k++) {
                // dans la spire, j = sigma * (-da/dt + Upq/Lz)   
                a[j + M * k][j + M * k] += coeff1; // -da/dt -> (aij(t-1) - aij(t))/dt
                a[j + M * k][M * N + i] = coeff2; // sigma*Upq/Lz
            }
        }
    }
    

    /* coefficients lies à l'integrale
     Ipq = sigma * integrale de (-da/dt + Upq/Lz)
         = sigma * somme de moyenne ponderee de (-da/dt + Upq/Lz)*dx*dy sur chaque carre du maillage (dans Spq)
         = ""                                "" ((a(t-1) - a(t))/dt + Upq/Lz)*dx*dy ""                      ""
         = ""                                "" ((a(t-1) - a(t))/dt)*dx*dy            + sigma*Upq/Lz * dx*dy*nombre de carres dans Spq
         = ""                                "" ((a(t-1) - a(t))/dt)*dx*dy            + sigma*Upq/Lz * surface de Spq

        ===>  sigma*sum(<-a(t)/dt>*dx*dy) + sigma*(Upq/Lz)*0.0002 - Ipq = sigma*sum(<-a(t-1)/dt>*dx*dy)
    */
    double n;
    coeff1 = sigma * dx * dy / dt;
    coeff2 = (sigma / Lz) * 0.0002;

    for (int i = 0; i < 4; i++) {
        int idx = M * N + i;
        int minx = ((int)(S[i][0] / dx)) - 1;
        int miny = ((int)(S[i][1] / dy)) - 1;
        int maxx = minx + ((int) (0.02 / dx));
        int maxy = miny + ((int)(0.01 / dx));
        // sigma* sum(<a(t) / dt>*dx * dy)
        for (int j = minx; j <= maxx; j++) {
            for (int k = miny; k <= maxy; k++) {
                /*     
                  n    = coeff lie à la moyenne ponderee dans l'equation lineaire
                  4*n  = nombre de carres du maillage de Spq dans lequel le point est inclu !

                              x---x---x---x         1---2---2---1
                              |   |   |   |         |   |   |   |
                  Spq  ===>   x---x---x---x   ===>  2---4---4---2 
                              |   |   |   |         |   |   |   |
                              x---x---x---x         1---2---2---1
                */
                if (j == minx || j == maxx || k == miny || k == maxy) {
                    if ((j == minx || j == maxx) && (k == miny || k == maxy)) n = 0.25; // 1/4 pour les coins
                    else n = 0.5; // 1/2 pour la frontiere (pas les coins)
                }
                else n = 1.0; // 1 pour les points interieurs (incidents à 4 carres du maillage)
                a[idx][j + M * k] = -n * coeff1; 
            }
        }
        a[idx][idx] = coeff2; //sigma* (Upq / Lz) * 0.0002
        a[idx][idx + 4] = -1.0; // -Ipq
    }

    // coefficients lies aux relations entre les spires
    a[M * N + 4][M * N] = 2.0; a[M * N + 4][M * N + 1] = 2.0;
    a[M * N + 5][M * N + 4] = 1.0; a[M * N + 5][M * N + 5] = -1.0;
    a[M * N + 6][M * N + 2] = 2.0; a[M * N + 6][M * N + 3] = 2.0; a[M * N + 6][M * N + 7] = -R; 
    a[M * N + 7][M * N + 6] = 1.0; a[M * N + 7][M * N + 7] = -1.0;
}

void* update_b(double* x, double* b, double t) {
    /* coefficients lies à j dans les spires
       dans la spire, j = sigma * (-da/dt + Upq/Lz)
       -da/dt -> (aij(t-1) - aij(t))/dt
    */
    double coeff = -sigma / dt;
    for (int i = 0; i < 4; i++) {
        int minx = ((int)(S[i][0] / dx)) - 1;
        int miny = ((int)(S[i][1] / dy)) - 1;
        int maxx = minx + ((int)(0.02 / dx));
        int maxy = miny + ((int)(0.01 / dx));
        for (int j = minx; j <= maxx; j++) {
            for (int k = miny; k <= maxy; k++) {
                b[j + k * M] = x[j + k * M] * coeff;
            }
        }
    }

    // coefficients lies à l'integrale
    // sigma*sum(<-a(t)/dt>*dx*dy) + sigma*(Upq/Lz)*0.0004 - Ipq = sigma*sum(<-a(t-1)/dt>*dx*dy)
    double n;
    coeff = sigma * dx * dy / dt;
    for (int i = 0; i < 4; i++) {
        int idx = M * N + i;
        int minx = ((int)(S[i][0] / dx)) - 1;
        int miny = ((int)(S[i][1] / dy)) - 1;
        int maxx = minx + ((int)(0.02 / dx));
        int maxy = miny + ((int)(0.01 / dx));
        b[idx] = 0; // reset à 0
        // sigma* sum(<a(t) / dt>*dx * dy)
        for (int j = minx; j <= maxx; j++) {
            for (int k = miny; k <= maxy; k++) {
                /*
                  n    = coeff lie à la moyenne pondérée dans l'equation lineaire
                  4*n  = nombre de carres du maillage de Spq dans lequel le point est inclu !

                              x---x---x---x         1---2---2---1
                              |   |   |   |         |   |   |   |
                  Spq  ===>   x---x---x---x   ===>  2---4---4---2
                              |   |   |   |         |   |   |   |
                              x---x---x---x         1---2---2---1
                */
                if (j == minx || j == maxx || k == miny || k == maxy) {
                    if ((j == minx || j == maxx) && (k == miny || k == maxy)) n = 0.25; // 1/4 pour les coins
                    else n = 0.5; // 1/2 pour la frontiere (pas les coins)
                }
                else n = 1.0; // 1 pour les points intérieurs (incidents à 4 carres du maillage)
                b[idx] -= x[j + M * k] * n;
            }
        }
        b[idx] *= coeff; // on multiplie par le coeff à le fin pour limiter les calculs inutiles!
    }

    // coefficient lie aux relations entre les spires
    b[M * N + 4] = -E1(t, f);
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        printf("Wrong number of arguments\n");
        return 0;
    }

    // initialisation de toutes les constantes utiles au bon fonctionnement du programme !
    L = 0.1;
    H = 0.1;
    Lz = 0.3;
    f = atof(argv[1]);
    M = atof(argv[2]);
    N = M;
    Mt = atof(argv[3]);
    S11[0] = 0.01 ; S11[1] =  0.055;
    S12[0] = 0.035; S12[1] =  0.055;
    S21[0] = 0.01 ; S21[1] =  0.035;
    S22[0] = 0.035; S22[1] =  0.035;
    dx = L / (M+1);
    dy = dx;
    dt = 1 / (Mt*f);
    mu0 = 4 * M_PI * 0.0000001;
    sigma = 5.97 * 10000000;
    R = 1.0;
    m = M*N+2*4;

    // creation des matrices et vecteurs utiles (remplis de 0 au depart)
    matrix* A = allocate_matrix(m, m);
    double* x = calloc(m, sizeof(double));
    double* b = calloc(m, sizeof(double));
    int* p = calloc(m, sizeof(int));
    double y[m];

    // initialisation du vecteur des pivots pour la factorisation LU
    for (int i = 0; i < m; i++) {
        p[i] = i;
    }
    
    FILE* output_file = fopen(argv[4], "w");
    // resolution du probleme

    // remplissage de la matrice A
    initialise_A(A);
    // decomposition LU avec pivots
    LU_decomposition_pivot(m, A, p);
    
    // resolution du systeme iterativement
    for (double t = 0; t <= 2/f; t += dt) {
        update_b(x, b, t); // mise a jour de b
        Forward_substitution_pivot(m, A, p, y, b); // forward substitution
        Back_substitution_pivot(m, A, p, x, y); // backward substitution
        fprintf(output_file, " %2.6e %2.6e %2.6e\n", t, E1(t, f) * x[M * M + 4], R * x[M * M + 7] * x[M * M + 7]); // impression des resultats dans le fichier    
    }
    
    // liberation des donnees et fermeture du fichier
    fclose(output_file);
    free(b);
    free(p);
    free(x);
    free_matrix(A);
    return 0 ;

}