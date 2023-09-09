#include <mpi/mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#define TEST_SEND
//#define ESTIMATE_SEND_TIME


int main(int argc, char* argv[]) {

    // Initializing MPI
    MPI_Init(&argc, &argv);

    // Getting the rank and communicator size
    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


// Test different modes of MPI_Send
#ifdef TEST_SEND
#define ARR_SIZE 3200000

    if(commsize != 2) {
        if(!rank)
            printf("Testing MPI_Send modes case requiers exactly 2 processes!\n");

        MPI_Finalize();
        return 0;
    }

    if(!rank)
        printf("Total sending size for MPI_Send tests in bytes: %ld\n", ARR_SIZE * sizeof(long long));

    // Create and fill a variable we will be sending and receiving
    long long* var = (long long*)calloc(ARR_SIZE, sizeof(long long));

    if(!rank)
        for(int i = 0; i < ARR_SIZE; i++)
            var[i] = i;

    // Ensure that both 0 and 1 processes are ready to perform given operations
    // and start at the same time to measure time correctly
    MPI_Barrier(MPI_COMM_WORLD);

    // Standard MPI_Send
    if(rank) {
        sleep(2);
        MPI_Recv(var, ARR_SIZE, MPI_LONG_LONG, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        double start = MPI_Wtime();
        MPI_Send(var, ARR_SIZE, MPI_LONG_LONG, rank + 1, 0, MPI_COMM_WORLD);
        double end   = MPI_Wtime();

        printf("Elapsed time due to MPI_Send:  %lg seconds\n", end - start);
    }

    // Ensure that both 0 and 1 processes are ready to perform given operations
    // and start at the same time to measure time correctly
    MPI_Barrier(MPI_COMM_WORLD);

    // Explicitly synchronous Send aka MPI_Ssend
    if(rank) {
        sleep(2);
        MPI_Recv(var, ARR_SIZE, MPI_LONG_LONG, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        double start = MPI_Wtime();
        MPI_Ssend(var, ARR_SIZE, MPI_LONG_LONG, rank + 1, 0, MPI_COMM_WORLD);
        double end   = MPI_Wtime();

        printf("Elapsed time due to MPI_Ssend: %lg seconds\n", end - start);
    }

    // Ensure that both 0 and 1 processes are ready to perform given operations
    // and start at the same time to measure time correctly
    MPI_Barrier(MPI_COMM_WORLD);

    // MPI_Rsend waits until matching Recv have posted
    if(rank) {
        sleep(2);
        MPI_Recv(var, ARR_SIZE, MPI_LONG_LONG, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        double start = MPI_Wtime();
        MPI_Rsend(var, ARR_SIZE, MPI_LONG_LONG, rank + 1, 0, MPI_COMM_WORLD);
        double end   = MPI_Wtime();

        printf("Elapsed time due to MPI_Rsend: %lg seconds\n", end - start);
    }

    // Ensure that both 0 and 1 processes are ready to perform given operations
    // and start at the same time to measure time correctly
    MPI_Barrier(MPI_COMM_WORLD);

    // Explicitly non-blocking Send aka MPI_Isend
    if(rank) {
        sleep(2);
        MPI_Recv(var, ARR_SIZE, MPI_LONG_LONG, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        // Have to pass request argument for MPI_Isend
        MPI_Request myRequest;

        double start = MPI_Wtime();
        MPI_Isend(var, ARR_SIZE, MPI_LONG_LONG, rank + 1, 0, MPI_COMM_WORLD, &myRequest);
        double end   = MPI_Wtime();

        printf("Elapsed time due to MPI_Isend: %lg seconds\n", end - start);
    }

    // Ensure that both 0 and 1 processes are ready to perform given operations
    // and start at the same time to measure time correctly
    MPI_Barrier(MPI_COMM_WORLD);

    // Bufferized Send aka MPI_Bsend
    if(rank) {
        sleep(2);
        MPI_Recv(var, ARR_SIZE, MPI_LONG_LONG, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        // This time we have to provide sufficient enough buffer for MPI_Bsend function
        // because default one is quite small for our purposes
        long long* tmpBuf = (long long*)calloc(ARR_SIZE, sizeof(long long));
        int tmpSize;

        MPI_Pack_size(ARR_SIZE, MPI_LONG_LONG, MPI_COMM_WORLD, &tmpSize);
        tmpSize += MPI_BSEND_OVERHEAD;
        MPI_Buffer_attach(tmpBuf, tmpSize);

        double start = MPI_Wtime();
        MPI_Bsend(var, ARR_SIZE, MPI_LONG_LONG, rank + 1, 0, MPI_COMM_WORLD);
        double end   = MPI_Wtime();

        printf("Elapsed time due to MPI_Bsend: %lg seconds\n", end - start);
        
        // Detaching our buffer and free its memory
        MPI_Buffer_detach(&tmpBuf, &tmpSize);
        free(tmpBuf);
    }    
    // Free allocated memory
    free(var);

#endif

// This section is for estimating time for P2P MPI_Send/MPI_Recv between processes
#ifdef ESTIMATE_SEND_TIME

    // Need to wait all processes because previous section can be executed as well
    // as this one at one time
    MPI_Barrier(MPI_COMM_WORLD);

    if(commsize != 2) {
        if(!rank)
            printf("Estimating send time case requiers exactly 2 processes!\n");

        MPI_Finalize();
        return 0;
    }

    /*
    To estimate time for sending/receiving, we:
    0) Create some array with large enough size
    
    1) Send some part of it, [array_size * sizeof(long long)] bytes, 
        as proc0 --> proc1 and receive it from proc1
    
    2) Evaluate time which was spent by these operations
    
    3) Repeat 1-2) [iters] times to get average time which was spent 
        by these operations with this fixed size
    
    4) Repeat 1-3) throughout [5, max_array_size) with [loop_step] step to obtain 
        average time to send message with different sizes
    */
    long int max_array_size = 12288;
    long long* sendArr = (long long*)calloc(max_array_size, sizeof(long long));

    double start = 0.0;
    double end   = 0.0;
    int iters = 100000;
    long int loop_step = 5;

    for(long int array_size = 5; array_size < max_array_size; array_size += loop_step) {
        double means = 0.0;

        for(int i = 0; i < iters; i++) {
        
            if(!rank) {
                start = MPI_Wtime();
                MPI_Send(sendArr, array_size, MPI_LONG_LONG, rank + 1, 0, MPI_COMM_WORLD);
            }
            else
                MPI_Recv(sendArr, array_size, MPI_LONG_LONG, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
            // We need to put some barriers to be sure that we measure whole data transfer
            MPI_Barrier(MPI_COMM_WORLD);

            if(!rank) {
                end   = MPI_Wtime();
                means += end - start;
            }

            // This one is for correctly calculating time within one [array_size]
            MPI_Barrier(MPI_COMM_WORLD);
        }

        // Output our estimations: current message byte size and average time to pass it
        if(!rank) {
            printf("%ld %.10lf\n", array_size * sizeof(long long), means / iters);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    // Free allocated memory
    free(sendArr);

#endif

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}