#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>


// Function to sum reciprocals
double Sum(long long start, long long end) {
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

    // Splitting default world communicator into two: one for proc with
    // rank = 0 and second for others
    MPI_Comm newComm;
    MPI_Comm_split(MPI_COMM_WORLD, !!rank, rank, &newComm);
    if(newComm == MPI_COMM_NULL) {
        printf("Process %d: could not split communicator\n", rank);
        MPI_Finalize();
        return 0;
    }

    // Now obtain new rank withing new communicators
    int newRank, newCommSize;
    MPI_Comm_rank(newComm, &newRank);
    MPI_Comm_size(newComm, &newCommSize);

    // Just print all information
    printf("Global rank:             %d\n"
           "Global communcator size: %d\n"
           "New rank:                %d\n"
           "New communcator size:    %d\n\n", rank, commsize, newRank, newCommSize);
    fflush(stdout);

    // Now we calculate [Sum 1/n] from n = 1 to n = N, but only inside commucator
    // which has multiple processes
    if(rank) {

        // Passing number which we should sum up to by argv[1]
        const long long N = atoll(argv[1]);

        // Distributing starts and ends of summing among processes
        const long long diff  = N / newCommSize;
        long long start = 1 + diff * newRank;
        long long end   = 1 + diff * (newRank + 1);

        if(N % newCommSize) {
            if(newRank < N % newCommSize) {
                start += newRank;
                end   += newRank + 1;
            }
            else {
                start += (N % newCommSize);
                end   += (N % newCommSize);
            }
        }

        // Calculate sum for each process
        double partSum = Sum(start, end);

        // Collecting all into 0-th process
        double totalSum;
        MPI_Reduce(&partSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, 0, newComm);
        if(!newRank) {
            printf("\nTotal sum S = %.16lg\n", totalSum);
            fflush(stdout);
        }
    }
    // And at the end we need to free our created communicator
    MPI_Comm_free(&newComm);

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}