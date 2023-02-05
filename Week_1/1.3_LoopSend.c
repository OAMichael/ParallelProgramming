#include <mpi/mpi.h>
#include <stdio.h>


int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int var = 0;

    if(rank)
        MPI_Recv(&var, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    printf("Process %d: current value = %d. Sending %d to process %d\n", rank, var, var + 1, (rank + 1) % commsize);

    var++;
    MPI_Send(&var, 1, MPI_INT, (rank + 1) % commsize, 0, MPI_COMM_WORLD);

    if(!rank) {
        MPI_Recv(&var, 1, MPI_INT, commsize - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process %d: current value = %d\n", rank, var);
    }

    MPI_Finalize();

    return 0;
}