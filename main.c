#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define m 1.0
#define alpha 1.0
#define MODE 0

void pochodne(double t, double *s, double *k, int N) {
    k[0] = 0;
    k[N] = 0;
    k[N+1] = 0;
    k[2*N+1] = 0;
    for(int i = 1; i < N-1; i++) {
        k[i] = s[N+i+1];
        k[N+1+i] = alpha/m*(s[i-1] -2*s[i] + s[i+1]);
    }
}

void rk4_vec(double t, double dt, int n, double *s, void (*f)(double , double *, double  *)) {
    #define M 1000
    static double k1[M], k2[M], k3[M], k4[M], w[M];
    int i;
    // Kopia tablicy
    for(i = 0; i < n; i++) {
        w[i] = s[i];
    }
    // get rk1
    f(t, w, k1);
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt/2*k1[i];
    }
    // get rk2
    f(t+dt/2, w, k2);
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt/2*k2[i];
    }
    // get rk3
    f(t+dt/2, w, k3);
    for(i = 0; i < n; i++) {
        w[i] = s[i] + dt*k3[i];
    }
    // get rk4
    f(t+dt, w, k4);
    // do s[]
    for(i = 0; i < n; i++) {
        s[i]= s[i] + dt/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
}

void changeSystem(double* s, double* s2) {
    s2[0] = s[0]*cos(s[1]);
    s2[1] = s[0]*sin(s[1]);
    s2[2] = s[2];
}

int main() {

    int N = 50;                     // ilość cząstek
    int n = 2*(N+1);                    // ilość zmiennych w układzie RRZ1
    double dt = 0.02;
    double nt = 5000;
    int tmax = (int)nt*dt;
    double delta = 0.1;
    double x_0 = 0;
    double x_max = N*delta;
    double t = 0;

    void (*f)(double, double *, double *, int);
    f = pochodne;

    double *s;
    s = (double *)malloc(n*sizeof(double));    // tablica  rozwiązań

    FILE *fp;
    FILE *fp2;

    switch(MODE){
        case 0:
            // plik z energiami
            fp = fopen("energie0.dat", "w");
            if(fp == NULL) {
                printf("Error: No file found\n");
                return EXIT_FAILURE;
            }

            // plik z polozeniami
            fp2 = fopen("polozenia.dat", "w");
            if(fp == NULL) {
                printf("Error: No file found\n");
                return EXIT_FAILURE;
            }

            // warunki początkowe, położenia
            for(int i = 0; i < N+1; i++) {
                double xi0 = delta*i;
                double sigma = 3*delta;
                s[i] = xi0 + delta/3.0f*exp( -pow((xi0 - x_max/2.0f), 2)/(2.0*pow(sigma, 2)));
            }

            // prędkości
            for(int i = N+2; i < 2*(N+1); i++) {
                s[i] = 0;
            }
            break;
    }
    
    double e_kin = 0;
    double e_pot = 0;
    double e = 0;

    // simulation loop
    for(size_t i = 1; i <= N; i++) {
        rk4_vec(t, dt, n, s, f);
        e_kin = 0;
        e_pot = 0;
        e = e_kin + e_pot;
        t=t+dt;
        fprintf(fp, "%f %f %f %f\n", t, e_kin, e_pot, e);
        for(int i = 0; i < N+1; i++) {
            fprintf(fp2, "%.2f", s[i]);
        }
        fprintf(fp2, "\n");
        
    }
    fclose(fp);
    fclose(fp2);
    return EXIT_SUCCESS;
}