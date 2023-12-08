#include <iostream>
#include <cmath>
#include <cstring>

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

    const int diff_J  = JSIZE / commsize;
    int start_J = diff_J * rank;
    int end_J   = diff_J * (rank + 1);

    if(JSIZE % commsize) {
        if(rank < JSIZE % commsize) {
            start_J += rank;
            end_J   += rank + 1;
        }
        else {
            start_J += (JSIZE % commsize);
            end_J   += (JSIZE % commsize);
        }
    }

    if (rank == commsize - 1) {
        end_J = JSIZE - 2;
    }

    for (int i = 0; i < ISIZE; ++i) {
        for (int j = 0; j < JSIZE; ++j) {
            a[i][j] = 10 * i + j;
        }
    }

    std::string outFilename = "task1_J_MPI_" + std::to_string(commsize) + ".dat";


    int localArraySize = end_J - start_J;
    double localArray[JSIZE];

    int* recvcnts = new int[commsize];
    int* displs = new int[commsize];
    for (int k = 0; k < commsize; ++k) {
        MPI_Gather(&localArraySize, 1, MPI_INT, recvcnts, 1, MPI_INT, k, MPI_COMM_WORLD);
    }

    displs[0] = 0;
    for(int i = 1; i < commsize; ++i) {
        displs[i] = displs[i - 1] + recvcnts[i - 1]; 
    }

    // Start time
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    if (commsize > 1) {
        for (int i = 3; i < ISIZE; ++i) {
            for (int j = start_J; j < end_J; ++j) {
                localArray[j] = std::sin(3 * a[i - 3][j + 2]);
            }

            MPI_Allgatherv(&localArray[start_J], localArraySize, MPI_DOUBLE, &a[i][0], recvcnts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        }
    }
    else {
        for (int i = 3; i < ISIZE; ++i) {
            for (int j = 0; j < JSIZE - 2; ++j) {
                a[i][j] = std::sin(3 * a[i - 3][j + 2]);
            }
        }
    }
    // End time
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();

    delete[] recvcnts;
    delete[] displs;

    if (!rank) {
        std::cout << "Number of executors: " << commsize << std::endl;
        std::cout << "Elapsed time: " << end - start << " seconds" << std::endl;

        serializeToFile(outFilename, a);
    }
    MPI_Finalize();
    return 0;
}
