#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static inline const int getNumDigits(const int num) {
    int tmp = num;
    unsigned numDigits = 0;
    while(tmp > 0) {
        numDigits++;
        tmp /= 10;
    }

    return numDigits;
}



int main(int argc, char* argv[]) {

    // Initializing MPI
    MPI_Init(&argc, &argv);

    // Getting the rank and communicator size
    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // 1 for '\0' and 1 for '\n'
    char* out = (char*)calloc(1 + 1 + getNumDigits(rank), sizeof(char));
    sprintf(out, "%d\n", rank);

    int out_chars = strlen(out);

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, "ordered_writing.txt", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    
    // Just for truncation
    MPI_File_set_size(fh, 0);
    
    MPI_File_write_at(fh, 2 * (commsize - 1 - rank), out, out_chars, MPI_CHAR, MPI_STATUS_IGNORE);


    free(out);
    MPI_File_close(&fh);
    
    // Finalizing MPI
    MPI_Finalize();

    return 0;
}