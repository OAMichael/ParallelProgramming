#include <iostream>
#include <cmath>

#include <mpi/mpi.h>

#include "Common.h"


int main(int argc, char **argv) {

    // Initializing MPI
    MPI_Init(&argc, &argv);

    // Getting the rank and communicator size
    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    My2DArray<double> a(ISIZE, JSIZE);

    const int diff_I  = ISIZE / commsize;
    int start_I = diff_I * rank;
    int end_I   = diff_I * (rank + 1);

    if(ISIZE % commsize) {
        if(rank < ISIZE % commsize) {
            start_I += rank;
            end_I   += rank + 1;
        }
        else {
            start_I += (ISIZE % commsize);
            end_I   += (ISIZE % commsize);
        }
    }

    for (int i = start_I; i < end_I; ++i) {
        for (int j = 0; j < JSIZE; ++j) {
            a[i][j] = 10 * i + j;
        }
    }

    std::string outFilename = "baseline_MPI_" + std::to_string(commsize) + ".dat";


    // Start time
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    if (commsize > 1) {
        
        for (int i = start_I; i < end_I; ++i) {
            for (int j = 0; j < JSIZE; ++j) {
                a[i][j] = std::sin(2 * a[i][j]);
            }
        }

        if (rank) {
            MPI_Request myRequest;

            int sendStart = start_I;            
            int sendSize = (end_I - start_I) * JSIZE;

            MPI_Isend(&sendStart, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &myRequest);
            MPI_Isend(&sendSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &myRequest);
            MPI_Isend(&a[start_I][0], sendSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &myRequest);
        }
        else {
            for (int k = 1; k < commsize; ++k) {
                int recvStart = 0;
                int recvSize = 0;

                MPI_Recv(&recvStart, 1, MPI_INT, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&recvSize, 1, MPI_INT, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&a[recvStart][0], recvSize, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    else {
        for (int i = 0; i < ISIZE; ++i) {
            for (int j = 0; j < JSIZE; ++j) {
                a[i][j] = std::sin(2 * a[i][j]);
            }
        }
    }
    // End time
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();


    if (!rank) {
        std::cout << "Number of executors: " << commsize << std::endl;
        std::cout << "Elapsed time: " << end - start << " seconds" << std::endl;

        serializeToFile(outFilename, a);
    }
    MPI_Finalize();
    return 0;
}
