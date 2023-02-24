#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmp/gmp.h"


typedef struct factorial {
    int n;
    mpz_t n_fact;
} factorial_t;




factorial_t CalculateMaxN(int N) {
    factorial_t tmp;

    mpz_init_set_ui(tmp.n_fact, 1);

    mpz_t inv_epsilon;
    mpz_init_set_ui(inv_epsilon, 10);
    mpz_pow_ui(inv_epsilon, inv_epsilon, N);

    int n = 0;
    while(mpz_cmp(inv_epsilon, tmp.n_fact) > 0) {
        n++;
        mpz_mul_ui(tmp.n_fact, tmp.n_fact, n);
    }

    mpz_clear(inv_epsilon);

    tmp.n = n;

    return tmp;
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


    // Passing number of accurate digits by argv[1]
    const int N = atoi(argv[1]);

    factorial_t MaxFact;
    if(rank == commsize - 1) {
        MaxFact = CalculateMaxN(N);

        for(int i = 0; i < commsize - 1; i++) {
            MPI_Send(&MaxFact.n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

    }
    else {
        MPI_Recv(&MaxFact.n, 1, MPI_INT, commsize - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    // Distributing starts and ends of summing among processes
    const int diff  = MaxFact.n / commsize;
    int start = 1 + diff * rank;
    int end   = 1 + diff * (rank + 1);

    if(MaxFact.n % commsize) {
        if(rank < MaxFact.n % commsize) {
            start += rank;
            end += rank + 1;
        }
        else {
            start += (MaxFact.n % commsize);
            end   += (MaxFact.n % commsize);
        }
    }


    mpz_t LocCurrFact;
    mpz_init_set_ui(LocCurrFact, 1);

    mpz_t LocSum;
    mpz_init_set_ui(LocSum, 1);

    for(int i = end - 1; i > start; i--) {
        mpz_mul_ui(LocCurrFact, LocCurrFact, i);
        mpz_add(LocSum, LocSum, LocCurrFact);
    }
    mpz_clear(LocCurrFact);


    mpz_t RankMaxFact;

    if(rank < commsize - 1)
        mpz_init_set_ui(RankMaxFact, 1);
    else
        mpz_init_set(RankMaxFact, MaxFact.n_fact);

    // For all except last processes we calculate RankMaxFact =  start * (start + 1) * ... * (end - 2) * (end - 1)
    // for future sending factorial to NEXT process
    if(rank < commsize - 1) {
        for(int i = start; i < end; i++) {
            mpz_mul_ui(RankMaxFact, RankMaxFact, i);
        }
    }


    char* FactStrRecv = NULL;

    if(rank && rank < commsize - 1) {
        // For all except first processes we receive largest factorial of PREVIOUS process
        MPI_Status status;
        int recvLength;
        MPI_Probe(rank - 1, 0, MPI_COMM_WORLD, &status);

        MPI_Get_count(&status, MPI_CHAR, &recvLength);

        FactStrRecv = (char*)calloc(recvLength, sizeof(char));

        MPI_Recv(FactStrRecv, recvLength, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // By these two lines RankFactFromStr is the largest factorial of PREVIOUS process 
        mpz_t RankFactFromStr;
        mpz_init_set_str(RankFactFromStr, FactStrRecv, 10);

        // If THIS process is not last, multiply largest factorial of PREVIOUS process by [start * (start + 1) * ... * (end - 2) * (end - 1)]
        // of THIS process to obtain largest factorial of THIS process
        mpz_mul(RankMaxFact, RankMaxFact, RankFactFromStr);

        mpz_clear(RankFactFromStr);
        free(FactStrRecv);
    }





    if(rank < commsize - 2) {
        char* FactStrSend = mpz_get_str(NULL, 10, RankMaxFact);
        MPI_Send(FactStrSend, strlen(FactStrSend) + 1, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD);
    }

    // By now all processes have value of their maximum factorial stored in RankMaxFact
    mpf_set_default_prec(64 + 8 * N);

    mpf_t LocSum_float;
    mpf_init(LocSum_float);
    mpf_set_z(LocSum_float, LocSum);

    mpf_t RankMaxFact_float;
    mpf_init(RankMaxFact_float);
    mpf_set_z(RankMaxFact_float, RankMaxFact);

    mpf_div(LocSum_float, LocSum_float, RankMaxFact_float);
    

    mpz_clear(RankMaxFact);
    mpz_clear(LocSum);
    if(rank == commsize - 1)
        mpz_clear(MaxFact.n_fact);

    mpf_clear(RankMaxFact_float);
    
    // Send sums from all processes to 0's process
    if(rank) {

        char *buf = (char*)calloc(N + 4, sizeof(char));

        char* formatStr = (char*)calloc(4 + strlen(argv[1]) + 1, sizeof(char));
        snprintf(formatStr, 4 + strlen(argv[1]) + 1, "%%.%dFf", N + 2);

        gmp_snprintf(buf, N + 4, formatStr, LocSum_float);
        free(formatStr);

        MPI_Send(buf, strlen(buf) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

        free(buf);
    }
    
    if(!rank) {
        mpf_add_ui(LocSum_float, LocSum_float, 1);

        // Now aggregate all sums from all processes
        for(int i = 1; i < commsize; i++) {

            MPI_Status status;
            int recvLength;
            MPI_Probe(i, 0, MPI_COMM_WORLD, &status);

            MPI_Get_count(&status, MPI_CHAR, &recvLength);

            char* strSum_i = (char*)calloc(recvLength, sizeof(char));

            MPI_Recv(strSum_i, recvLength, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            
            mpf_t Sum_i;
            mpf_init_set_str(Sum_i, strSum_i, 10);
            free(strSum_i);

            mpf_add(LocSum_float, LocSum_float, Sum_i);
            mpf_clear(Sum_i);
            
        }

        
        char* formatStr = (char*)calloc(6 + strlen(argv[1]) + 1, sizeof(char));

        snprintf(formatStr, 6 + strlen(argv[1]), "%%.%dFf\n", N);
        gmp_printf(formatStr, LocSum_float);

        free(formatStr);
        mpf_clear(LocSum_float);
        
    }

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}