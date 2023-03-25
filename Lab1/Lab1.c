#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


//#define ONE_PROC_FAST

const double a = 2.0;

const double T = 1;
const double X = 1;

const int K = 500 + 1;
const int M = 2000 + 1;

const double tau = T / (K - 1);
const double h   = X / (M - 1);



double phi(double x) {
    return cos(3.14159 * x);
}


double psi(double t) {
    return exp(-t);
}


double f(double x, double t) {
    return x + t;
}




/*
double phi(double x) {
    return cos(3.14159 * x);
}


double psi(double t) {
    return exp(-t);
}


double f(double x, double t) {
    return 2 + exp(-x) * sin(t) + exp(-2 * x) * sin(2 * t) + exp(-4 * x) * sin(4 * t) +
               cos(0.0135 * x) / atan(1.112 * (1.0 + t)) + sinh(x * cosh(x * t)) + 1 / (1 - t * x / 2.0);
}

*/


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


    int iterations = 10000;
    double avg = 0;

#ifdef ONE_PROC_FAST

    for(int iter = 0; iter < iterations; iter++) {

        double start_one = MPI_Wtime();
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
                u[k][m] = 2 * tau * h / (tau + h) * f((m + 0.5) * h, (k - 0.5) * tau) + u[k - 1][m - 1] + (u[k][m - 1] - u[k - 1][m]) * (tau - h) / (tau + h);
            }
        }

        printf("t x u\n");

        for(int k = 0; k < K; k++) {
            //for(int m = 0; m < M; m++) {
                //printf("%.6lf %.6lf %.6lf\n", k * tau, m * h, u[k][m]);
            //}
            free(u[k]);
        }
        free(u);
        


        double end_one = MPI_Wtime();

        avg += end_one - start_one;

    }
    printf("Average time for one proc: %lf\n", avg / iterations);

#else

    long int tmp = 0;

    double startTime = 0;
    double endTime   = 0;

    //double* tmpBuf = (double*)calloc(2 * iterations * (K - 1), sizeof(double));
    //int tmpSize;

    //MPI_Pack_size(2 * iterations * (K - 1), MPI_DOUBLE, MPI_COMM_WORLD, &tmpSize);
    //tmpSize += MPI_BSEND_OVERHEAD * 2 * iterations * (K - 1);
    //MPI_Buffer_attach(tmpBuf, tmpSize);

    for(int iter = 0; iter < iterations; iter++) {

        if(!rank)
            startTime = MPI_Wtime();

        // Distributing starts and ends along M among processes
        const long long diff_M  = M / commsize;
        long long start_M = diff_M * rank;
        long long end_M   = diff_M * (rank + 1);

        if(M % commsize && rank == commsize - 1) {
            if(rank < M % commsize) {
                start_M += rank;
                end_M   += rank + 1;
            }
            else {
                start_M += (M % commsize);
                end_M   += (M % commsize);
            }
        }


        double** u = (double**)calloc(K, sizeof(double*));

        for(int i = 0; i < K; i++) {
            u[i] = (double*)calloc(M, sizeof(double));
        }


        //if(!rank)    
            //for(int k = 0; k < K; k++)
                //u[k][0] = psi(k * tau);

        //for(int m = start_M - !!rank; m < end_M; m++)
            //u[0][m] = phi(m * h);


        //const double c_1 =  2 * tau * h  / (a * tau + h);
        //const double c_2 = (a * tau - h) / (a * tau + h);
        //MPI_Request myRequest;


        if(rank)
            for(int k = 1; k < K; k++) {
                //MPI_Recv(&u[k][start_M - 1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for(int m = start_M; m < end_M; m++) {
                    //u[k][m] = c_1 * f((m + 0.5) * h, (k - 0.5) * tau) + u[k - 1][m - 1] + (u[k][m - 1] - u[k - 1][m]) * c_2;
                    tmp += m % 2;
                }

                //if(rank < commsize - 1) {
                    //MPI_Bsend(&u[k][end_M - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                //}
            }
        else 
            for(int k = 1; k < K; k++) {
                for(int m = start_M + 1; m < end_M; m++) {
                    //u[k][m] = c_1 * f((m + 0.5) * h, (k - 0.5) * tau) + u[k - 1][m - 1] + (u[k][m - 1] - u[k - 1][m]) * c_2;
                    tmp += m % 2;
                }

                //if(commsize > 1) {
                    //MPI_Bsend(&u[k][end_M - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                //}
            }

        //printf("%lf\n", tmp);

        //MPI_Barrier(MPI_COMM_WORLD);
        //if(!rank)
            //printf("t x u\n");

        //int sent = 0;
        //if(rank)
            //MPI_Recv(&sent, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for(int k = 0; k < K; k++) {
            //for(int m = start_M; m < end_M; m++) {
                //printf("%.6lf %.6lf %.6lf\n", k * tau, m * h, u[k][m]);
            //}
            free(u[k]);
        }
        free(u);

        //if(rank < commsize - 1) {
            //MPI_Request myRequest;
            //MPI_Isend(&sent, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &myRequest);
        //}

        //break;
        MPI_Barrier(MPI_COMM_WORLD);

        if(!rank) {
            endTime = MPI_Wtime();
            avg += endTime - startTime;
        }
    }
    if(!rank)
        printf("Average time for %d processes: %lf\n", commsize, avg / iterations);
    printf("%ld\n", tmp);

    //MPI_Buffer_detach(&tmpBuf, &tmpSize);
    //free(tmpBuf);
    



    /*
    double startTime = 0;
    double endTime   = 0;
    for(int iter = 0; iter < iterations; iter++) {

        if(!rank)
            startTime = MPI_Wtime();

        // Distributing starts and ends along M among processes
        const long long diff_M  = M / commsize;
        long long start_M = diff_M * rank;
        long long end_M   = diff_M * (rank + 1);

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

        int total = K * M;

        double* u = (double*)calloc(total, sizeof(double));


        if(!rank)    
            for(int k = 0; k < K; k++)
                u[k * M] = psi(k * tau);

        for(int m = start_M; m < end_M; m++)
            u[m] = phi(m * h);


        double c_1 =  2 * tau * h  / (a * tau + h);
        double c_2 = (a * tau - h) / (a * tau + h);

        //double startMainLoop = MPI_Wtime();
        for(int k = 1; k < K; k++)
        {
            if(rank) {
                MPI_Recv(&u[k * M + start_M - 1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for(int m = start_M; m < end_M; m++) {
                    u[k * M + m] = c_1 * f((m + 0.5) * h, (k - 0.5) * tau) + u[(k - 1) * M + m - 1] + (u[k * M + m - 1] - u[(k - 1) * M + m]) * c_2;
                }

                if(rank < commsize - 1) {
                    MPI_Bsend(&u[k * M + end_M - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                }
                continue;
            }

            for(int m = start_M + 1; m < end_M; m++) {
                    u[k * M + m] = c_1 * f((m + 0.5) * h, (k - 0.5) * tau) + u[(k - 1) * M + m - 1] + (u[k * M + m - 1] - u[(k - 1) * M + m]) * c_2;
            }

            if(commsize > 1) {
                MPI_Bsend(&u[k * M + end_M - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            }
        }

        //double endMainLoop = MPI_Wtime();
        //printf("Process %d: elapsed time %lf\n", rank, endMainLoop - startMainLoop);

        //if(!rank)
            //printf("t x u\n");

        int sent = 0;
        //if(rank)
            //MPI_Recv(&sent, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for(int k = 0; k < total; k++) {
            //for(int m = start_M; m < end_M; m++) {
                sent = rank;
                //printf("%.6lf %.6lf %.6lf\n", k * tau, m * h, u[k][m]);
            //}
            //free(u[k]);
        }
        free(u);

        //if(rank < commsize - 1) {
            //MPI_Request myRequest;
            //MPI_Isend(&sent, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &myRequest);
        //}

        MPI_Barrier(MPI_COMM_WORLD);

        if(!rank) {
            endTime = MPI_Wtime();
            avg += endTime - startTime;
        }
    }
    if(!rank)
        printf("Average time for %d processes: %lf\n", commsize, avg / iterations);
    */
#endif

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}