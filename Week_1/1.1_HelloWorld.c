#include <mpi/mpi.h>
#include <stdio.h>


int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hello World! Communicator size = %d and my rank = %d\n", commsize, rank);

    MPI_Finalize();

    return 0;
}