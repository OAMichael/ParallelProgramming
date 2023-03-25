#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>


// Function to sum reciprocals
static inline double Sum(const long long start, const long long end) {
    double partSum = 0.0;

    for(long long i = start; i < end; i++) {
        partSum += 1.0 / (double)i;
    }

    return partSum;
}



int main(int argc, char* argv[]) {

    // Initializing MPI
    MPI_Init(&argc, &argv);

    // Getting the rank and communicator size
    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // In case of wrong input
    if(argc < 2) {
        if(!rank)
            printf("Usage: %s [N]\n", argv[0]);

        MPI_Finalize();
        return 0;
    }

    // Passing number which we should sum up to by argv[1]
    const long long N = atoll(argv[1]);

    // Distributing starts and ends of summing among processes
    const long long diff  = N / commsize;
    long long start = 1 + diff * rank;
    long long end   = 1 + diff * (rank + 1);

    if(N % commsize) {
        if(rank < N % commsize) {
            start += rank;
            end += rank + 1;
        }
        else {
            start += (N % commsize);
            end   += (N % commsize);
        }
    }

    // Calculate sum for each process
    double partSum = Sum(start, end);
    printf("Process %d: summing interval = [%10lld, %10lld)  ==>  Sum = %.16lg\n", rank, start, end, partSum);

    /*
    int MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win)
    */

    MPI_Win win;
    MPI_Win_create(&partSum, (MPI_Aint)sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);


    /* 
    int MPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype, 
                       int target_rank, MPI_Aint target_disp, int target_count, 
                       MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
    */

    if(rank) {   
        MPI_Accumulate(&partSum, 1, MPI_DOUBLE, 0, (MPI_Aint)0, 1, MPI_DOUBLE, MPI_SUM, win);
    }

    MPI_Win_fence(0, win);

    if(!rank)
        printf("Total sum S = %.16lg\n", partSum);


    MPI_Win_free(&win);

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}