/*
              / du/dt + a * du/dx = f(x, t),     0 <= t <= T, 0 <= x <= X
    We solve <| u(x, 0) = phi(x),                0 <= x <= X
              \ u(0, t) = psi(t),                0 <= t <= T



    To solve the PDE let us introduce a discrete mesh:

        t = k * tau, 0 <= k <= K
        x = m * h,   0 <= m <= M

        K * tau = T
        M * h   = X

    
    We will be using 2D explicit finite difference method (FDM):

        (u[k+1][m-1] - u[k][m-1] + u[k+1][m] - u[k][m]) / (2*tau) + (u[k+1][m] - u[k+1][m-1] + u[k][m] - u[k][m-1]) / (2*h) = f[k+1/2][m+1/2]

    where k and m correspond to time and space steps respectively. 
    Or, if we change k to start from 1:

        (u[k][m-1] - u[k-1][m-1] + u[k][m] - u[k-1][m]) / (2*tau) + (u[k][m] - u[k][m-1] + u[k-1][m] - u[k-1][m-1]) / (2*h) = f[k-1/2][m+1/2]


    It can be visualized so (assuming t axis goes up and x right):
    
   (k, m-1)  _         _  (k, m)
            |_| ----> |_|
                   __
                    /\ ^
                   /   |
                  /    |
                 /     |
                /      |    
               /       |
             _         _
            |_|       |_|
  (k-1, m-1)              (k-1, m)

    Thus, we start from left bottom corner of the mesh and every step evaluate (k, m) node from (k-1, m-1), (k-1, m) and (k, m-1) nodes.

    Sequential algorithm is straight forward: we will go row by row until work is done. Every step requires 
    knowledge of (k-1, m-1), (k-1, m) and (k, m-1) nodes. But since we go row by row, we already know these values
    for each step.
    
      t |_____________________________________
        |     |     |     |     |     |     |_
        |  >  |  >  |  >  |  >  |  >  |  >  |_
        |_____|_____|_____|_____|_____|_____|_
        |     |     |     |     |     |     |_
    K   |  >  |  >  |  >  |  >  |  >  |  >  |_
        |_____|_____|_____|_____|_____|_____|_
        |     |     |     |     |     |     |_
        |  >  |  >  |  >  |  >  |  >  |  >  |_
        |_____|_____|_____|_____|_____|_____|__
       0                                          x
                          M



    Parallel algorithm is quite interesting. We divide space line into several regions, one for each process. Within its region
    each process performs sequential algorithm, but as far as it on the last node on, it sends this value to next process. Only 
    after this operation it starts new row. This enables next process to start its corresponding row. The algorithm goes on 
    until work is done. Whole mesh looks like this:
    
     
      t |
        |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_
        |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_
        |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_
        |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_
     K  |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_
        |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_
        |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_
        |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_
        |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|_| |_|_|_|__
       0  proc0     proc1     proc2     proc3     ...    x
                                M


    And somewhere in the middle of completion it looks like this:

      t |
        |_____|_____|_____| |_____|_____|_____| |_____|_____|_____| |_____|_____|_____|___
        |     |     |     | |     |     |     | |     |     |     | |     |     |     |
        |     |     |     | |     |     |     | |     |     |     | |     |     |     |
        |_____|_____|_____| |_____|_____|_____| |_____|_____|_____| |_____|_____|_____|___
        |     |     |     | |     |     |     | |     |     |     | |     |     |     |
        |  >  | ... |     | |     |     |     | |     |     |     | |     |     |     |
        |_____|_____|_____| |_____|_____|_____| |_____|_____|_____| |_____|_____|_____|___
        |     |     |     | |     |     |     | |     |     |     | |     |     |     |
        |  >  |  >  |  >  | |  >  |  >  | ... | |     |     |     | |     |     |     |
        |_____|_____|_____| |_____|_____|_____| |_____|_____|_____| |_____|_____|_____|___
        |     |     |     | |     |     |     | |     |     |     | |     |     |     |
        |  >  |  >  |  >  | |  >  |  >  |  >  | |  >  |  >  |  >  | |     |     |     |
        |_____|_____|_____| |_____|_____|_____| |_____|_____|_____| |_____|_____|_____|___
        |     |     |     | |     |     |     | |     |     |     | |     |     |     |
        |  >  |  >  |  >  | |  >  |  >  |  >  | |  >  |  >  |  >  | |  >  |  >  | ... |
        |_____|_____|_____| |_____|_____|_____| |_____|_____|_____| |_____|_____|_____|___
        |     |     |     | |     |     |     | |     |     |     | |     |     |     |
        |  >  |  >  |  >  | |  >  |  >  |  >  | |  >  |  >  |  >  | |  >  |  >  |  >  |
        |_____|_____|_____| |_____|_____|_____| |_____|_____|_____| |_____|_____|_____|_____
       0                                                                                    x
               proc0               proc1               proc2               proc3        ...      

    
    Parallel algorithm will be implemented using MPI (Message Passing Interface).

    Note that in case of only 1 process the parallel algorithm degenerates to sequential 
    one without even overhead costs. So, it becomes purely sequential algorithm.
*/


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


// Define some problem specific constants
#define K 401
#define M 401


const double a = 2.0;

const double T = 1;
const double X = 1;

const double tau = T / (K - 1);
const double h   = X / (M - 1);


// Boundary and initial fuctions
double phi(double x) {
    return cos(3.14159 * x);
}


double psi(double t) {
    return exp(-t);
}


// RHS function
double f(double x, double t) {
    return x + t;
}



int main(int argc, char* argv[]) {

    // Initializing MPI
    MPI_Init(&argc, &argv);

    // Getting the rank and communicator size
    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // To speed up a little process interactions, we will be using non-blocking Bsend function 
    // Therefore, we need to provide sufficient buffer for it
    int Bsend_Size;

    MPI_Pack_size(2 * (K - 1), MPI_DOUBLE, MPI_COMM_WORLD, &Bsend_Size);
    Bsend_Size += MPI_BSEND_OVERHEAD * 2 * (K - 1);
    double* Bsend_Buf = (double*)malloc(Bsend_Size);

    MPI_Buffer_attach(Bsend_Buf, Bsend_Size);

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

    // This is a total width of process region in terms of nodes
    long long total_M = end_M - start_M;

    // For memory keeping we will be using only previous and current rows at one moment
    double* prev_row = (double*)calloc(total_M, sizeof(double));
    double* curr_row = (double*)calloc(total_M, sizeof(double));

    // This is done to properly plot a 3D graph of solution later
    if(!rank)
        printf("t x u\n");

    MPI_Barrier(MPI_COMM_WORLD);

    // For each process initialize 0-th row with corresponding phi(x) value
    for(int m = 0; m < total_M; m++) {
        prev_row[m] = phi((start_M + m) * h);
        printf("%.6lf %.6lf %.6lf\n", 0.0f, (start_M + m) * h, prev_row[m]);
    }

    // Calculate those coefficients once
    const double c_1 =  2 * tau * h  / (a * tau + h);
    const double c_2 = (a * tau - h) / (a * tau + h);

    // Temporal pointer is used only for row pointers swap
    double* p_tmp = NULL;

    // Main part of the algorithm
    if(rank) {
        // For additional processes
        double prev_col_curr_row = 0;
        double prev_col_prev_row = phi((start_M - 1) * h);

        for(int k = 1; k < K; k++) {
            // Receive value of (current row, previous column) as for begining new row 
            MPI_Recv(&prev_col_curr_row, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Perform all compuations according to the FDM we proposed at the begining
            curr_row[0] = c_1 * f((start_M + 0.5) * h, (k - 0.5) * tau) + prev_col_prev_row + (prev_col_curr_row - prev_row[0]) * c_2;
            printf("%.6lf %.6lf %.6lf\n", k * tau, start_M * h, curr_row[0]);

            for(int m = 1; m < total_M; m++) {
                curr_row[m] = c_1 * f((start_M + m + 0.5) * h, (k - 0.5) * tau) + prev_row[m - 1] + (curr_row[m - 1] - prev_row[m]) * c_2;
                printf("%.6lf %.6lf %.6lf\n", k * tau, (start_M + m) * h, curr_row[m]);
            }

            // Swap rows' pointers
            p_tmp = curr_row;
            curr_row = prev_row;
            prev_row = p_tmp;

            // Since for new compuation we need not only value of (current row, previous column) but
            // (previous row, previous column) as well, we store 'current' as 'previous' for next iteration
            prev_col_prev_row = prev_col_curr_row;

            // Send last value of row to next process to continue this chain
            if(rank < commsize - 1) {
                MPI_Bsend(&prev_row[total_M - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            }
        }
    }
    else {
        // For 0-th, main, process
        for(int k = 1; k < K; k++) {
            // Perform all compuations according to the FDM we proposed at the begining
            curr_row[0] = psi(k * tau);
            printf("%.6lf %.6lf %.6lf\n", k * tau, 0.0, curr_row[0]);

            for(int m = 1; m < total_M; m++) {
                curr_row[m] = c_1 * f((m + 0.5) * h, (k - 0.5) * tau) + prev_row[m - 1] + (curr_row[m - 1] - prev_row[m]) * c_2;
                printf("%.6lf %.6lf %.6lf\n", k * tau, m * h, curr_row[m]);
            }

            // Swap rows' pointers
            p_tmp = curr_row;
            curr_row = prev_row;
            prev_row = p_tmp;

            // Send last value of row to next process to continue this chain
            if(commsize > 1) {
                MPI_Bsend(&prev_row[total_M - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            }
        }
    }

    // Free used memory and detach buffer used by Bsend function
    free(curr_row);
    free(prev_row);

    free(Bsend_Buf);
    MPI_Buffer_detach(&Bsend_Buf, &Bsend_Size);


    // Finalizing MPI
    MPI_Finalize();

    return 0;
}