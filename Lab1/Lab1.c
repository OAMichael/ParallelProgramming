#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


//#define ONE_PROC_FAST

const double a = 2.0;

const double T = 1;
const double X = 1;

const int K = 1200 + 1;
const int M = 2000 + 1;

const double tau = T / (K - 1);
const double h   = X / (M - 1);


double phi(double x) {
    return cos(3.1415926536 * x);
}


double psi(double t) {
    return exp(-t);
}


double f(double x, double t) {
    return x + t;
}



/*
  t |
    |_|_|_|_|_|_|_|_
    |_|_|_|_|_|_|_|_
 K  |_|_|_|_|_|_|_|_
    |_|_|_|_|_|_|_|_
    |_|_|_|_|_|_|_|_
    |_|_|_|_|_|_|_|__
    0                x
            M
*/

int main(int argc, char* argv[]) {

    // Initializing MPI
    MPI_Init(&argc, &argv);

    // Getting the rank and communicator size
    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


#ifdef ONE_PROC_FAST

    double** u = (double**)calloc(K, sizeof(double*));

    for(int i = 0; i < K; i++) {
        u[i] = (double*)calloc(M, sizeof(double));
    }

    for(int k = 0; k < K; k++)
        u[k][0] = psi(k * tau);

    for(int m = 0; m < M; m++)
        u[0][m] = phi(m * h);


    for(int k = 1; k < K; k++) {
        for(int m = 1; m < M; m++) {
            u[k][m] = 2 * tau * h / (a * tau + h) * f((m + 0.5) * h, (k - 0.5) * tau) + u[k - 1][m - 1] + (u[k][m - 1] - u[k - 1][m]) * (a * tau - h) / (a * tau + h);
        }
    }

    printf("t x u\n");

    for(int k = 0; k < K; k++) {
        for(int m = 0; m < M; m++) {
            printf("%.6lf %.6lf %.6lf\n", k * tau, m * h, u[k][m]);
        }
        free(u[k]);
    }
    free(u);

#else

    // Distributing starts and ends along M among processes
    int diff_M  = M / commsize;
    int start_M = diff_M * rank;
    int end_M   = diff_M * (rank + 1);

    if(M % commsize) {
        if(rank < M % commsize) {
            start_M += rank;
            end_M   += rank + 1;
        }
        else {
            start_M += (M % commsize);
            end_M   += (M % commsize);
        }
    }

    diff_M = end_M - start_M;


    double** u = (double**)calloc(K, sizeof(double*));
    for(int i = 0; i < K; i++) {
        u[i] = (double*)calloc(diff_M, sizeof(double));
    }
    
    if(!rank)
        for(int k = 0; k < K; k++)
            u[k][0] = psi(k * tau);


    for(int m = 0; m < diff_M; m++)
        u[0][m] = phi((m + start_M) * h);



    const double a_1 = 2 * tau * h / (a * tau + h);
    const double a_2 = (a * tau - h) / (a * tau + h);

    double new_u  = 0;
    double prev_u = phi((start_M - 1) * h);

    for(int k = 1; k < K; k++)
    {
        if(rank) {
            MPI_Recv(&new_u, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            u[k][0] = a_1 * f((start_M + 0.5) * h, (k - 0.5) * tau) + prev_u + (new_u - u[k - 1][0]) * a_2;


            for(int m = 1; m < diff_M; m++) {
                u[k][m] = a_1 * f((m + start_M + 0.5) * h, (k - 0.5) * tau) + u[k - 1][m - 1] + (u[k][m - 1] - u[k - 1][m]) * a_2;
            }

            if(rank < commsize - 1)
                MPI_Send(&u[k][diff_M - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

            prev_u = new_u;
        }
        else {
            for(int m = 1; m < diff_M; m++) {
                u[k][m] = a_1 * f((m + start_M + 0.5) * h, (k - 0.5) * tau) + u[k - 1][m - 1] + (u[k][m - 1] - u[k - 1][m]) * a_2;
            }

            if(commsize > 1)
                MPI_Send(&u[k][diff_M - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }

    }


    if(!rank)
        printf("t x u\n");

    int sent = 0;
    if(rank)
        MPI_Recv(&sent, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for(int k = 0; k < K; k++) {
        for(int m = 0; m < diff_M; m++) {

            printf("%.6lf %.6lf %.6lf\n", k * tau, (start_M + m) * h, u[k][m]);
        }
        free(u[k]);
    }
    free(u);

    if(rank < commsize - 1)
        MPI_Send(&sent, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);

    
#endif

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}