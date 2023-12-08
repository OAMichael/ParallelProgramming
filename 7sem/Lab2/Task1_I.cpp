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

    for (int i = 0; i < ISIZE; ++i) {
        for (int j = 0; j < JSIZE; ++j) {
            a[i][j] = 10 * i + j;
        }
    }


    std::string outFilename = "task1_I_MPI_" + std::to_string(commsize) + ".dat";


    // Start time
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    if (commsize > 2) {
        const int end_I = ISIZE - (ISIZE - 3) % commsize;
        for (int i = 3 + rank; i < end_I; i += commsize) {

            if (commsize != 3 && i >= 6) {
                MPI_Recv(&a[i - 3][0], JSIZE, MPI_DOUBLE, (rank - 3 + commsize) % commsize, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            for (int j = 0; j < JSIZE - 2; ++j) {
                a[i][j] = std::sin(3 * a[i - 3][j + 2]);
            }
    
            MPI_Request myRequest;
            if (commsize != 3) {
                MPI_Isend(&a[i][0], JSIZE, MPI_DOUBLE, (rank + 3) % commsize, 0, MPI_COMM_WORLD, &myRequest);
            }

            if (rank) {
                MPI_Isend(&a[i][0], JSIZE, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &myRequest);
            }
            else {
                for (int k = 1; k < commsize; ++k) {
                    MPI_Recv(&a[i + k][0], JSIZE, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
                }
            }
        }

        // Handle some corner cases
        if (!rank) {
            for (int i = end_I; i < ISIZE; ++i) {
                for (int j = 0; j < JSIZE - 2; ++j) {
                    a[i][j] = std::sin(3 * a[i - 3][j + 2]);
                }
            }
        }
    }
    else if (commsize == 2) {
        if (!rank) {
            // 0-th rank. Calculate 3, 6, 9, 12, ... th rows
            const int end_J = ISIZE - ISIZE % 3;
            for (int i = 3; i < ISIZE; i += 3) {
                for (int j = 0; j < JSIZE - 2; ++j) {
                    a[i][j] = std::sin(3 * a[i - 3][j + 2]);
                }
            }

            for (int i = 4; i < end_J; i += 3) {
                MPI_Recv(&a[i][0], 2 * JSIZE, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
            }

            // Handle some corner cases
            const int remain = ISIZE % 3;
            for (int i = 0; i < remain; ++i) {
                for (int j = 0; j < JSIZE - 2; ++j) {
                    a[ISIZE - 1 - i][j] = std::sin(3 * a[ISIZE - 4 - i][j + 2]);
                }
            }         
        }
        else {
            // 1-st rank. Calculate 4, 5, 7, 8, 10, 11, ... th rows
            const int end_J = ISIZE - ISIZE % 3;
            for (int i = 4; i < end_J; i += 3) {
                for (int j = 0; j < JSIZE - 2; ++j) {
                    a[i][j] = std::sin(3 * a[i - 3][j + 2]);
                    a[i + 1][j] = std::sin(3 * a[i - 2][j + 2]);
                }
                MPI_Request myRequest;
                MPI_Isend(&a[i][0], 2 * JSIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &myRequest);
            }
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


    if (!rank) {
        std::cout << "Number of executors: " << commsize << std::endl;
        std::cout << "Elapsed time: " << end - start << " seconds" << std::endl;

        serializeToFile(outFilename, a);
    }
    MPI_Finalize();
    return 0;
}
