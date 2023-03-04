#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "gmp/gmp.h"


static inline int CalculateMaxN(const int N) {

    double x_curr = 3.0;
    double x_prev = x_curr;

    // x_{n + 1} = x_n - f(x_n)/f'(x_n),
    // where f(x) = x*ln(x) - x - N*ln(10)
    do {
        x_prev = x_curr;
        x_curr = (x_curr + N * log(10)) / log(x_curr);
    } while(fabs(x_curr - x_prev) > 1.0);

    return (int)ceil(x_curr);
}


static inline void mpz_set_ull(mpz_t n, int_fast64_t ull)
{
    mpz_set_ui(n, (unsigned int)(ull >> 32));
    mpz_mul_2exp(n, n, 32);
    mpz_add_ui(n, n, (unsigned int)ull);
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


    // Calculate to which term we must calculate for given accuracy and send to all processes this information
    // Still valid if we have only 1 process
    int MaxFact;
    if(rank == commsize - 1)
        MaxFact = CalculateMaxN(N);

    MPI_Bcast(&MaxFact, 1, MPI_INT, commsize - 1, MPI_COMM_WORLD);

    // Distributing starts and ends of summing among processes
    const int diff  = MaxFact / commsize;
    int start = 1 + diff * rank;
    int end   = 1 + diff * (rank + 1);

    if(MaxFact % commsize) {
        if(rank < MaxFact % commsize) {
            start += rank;
            end += rank + 1;
        }
        else {
            start += (MaxFact % commsize);
            end   += (MaxFact % commsize);
        }
    }
    mpz_t a_mpz;
    mpz_init_set_ui(a_mpz, end - 1);

    mpz_t S_mpz;
    mpz_init_set_ui(S_mpz, end);

    // Main algorithm
    mpz_t LocCurrFact;
    mpz_init_set_ui(LocCurrFact, 1);

    mpz_t LocSum;
    mpz_init_set_ui(LocSum, 0);


    for(int i = end - 2; i > start; i--) {
        if(i % 8192 == 0) {
            mpz_addmul(LocSum, LocCurrFact, S_mpz);
            mpz_mul(LocCurrFact, LocCurrFact, a_mpz);

            mpz_set_ull(a_mpz, i);
            mpz_set(S_mpz, a_mpz);
        }
        else {
            mpz_mul_ui(a_mpz, a_mpz, i);
            mpz_add(S_mpz, S_mpz, a_mpz);
        }
    }

    mpz_addmul(LocSum, LocCurrFact, S_mpz);
    mpz_mul(LocCurrFact, LocCurrFact, a_mpz);
    mpz_clears(a_mpz, S_mpz, NULL);
    

    // For all we calculate LocCurrFact = start * (start + 1) * ... * (end - 2) * (end - 1)
    mpz_mul_ui(LocCurrFact, LocCurrFact, start);


    if(rank) {
        // For all except first processes we receive largest factorial of PREVIOUS process
        MPI_Status status;
        int recvLength;

        // Use MPI_Probe and MPI_Get_count to obtain length of passed number
        MPI_Probe(rank - 1, 0, MPI_COMM_WORLD, &status);

        MPI_Get_count(&status, MPI_CHAR, &recvLength);

        char* FactStrRecv = (char*)calloc(recvLength, sizeof(char));

        MPI_Recv(FactStrRecv, recvLength, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // By these two lines RankFactFromStr is the largest factorial of PREVIOUS process 
        mpz_t RankFactFromStr;
        mpz_init_set_str(RankFactFromStr, FactStrRecv, 32);

        // If THIS process is not last, multiply largest factorial of PREVIOUS process by [start * (start + 1) * ... * (end - 2) * (end - 1)]
        // of THIS process to obtain largest factorial of THIS process
        mpz_mul(LocCurrFact, LocCurrFact, RankFactFromStr);

        mpz_clear(RankFactFromStr);
        free(FactStrRecv);
    }

    // Send NEXT process largest factorial of THIS process
    // Still valid if we have only 1 process
    if(rank < commsize - 1) {
        char* FactStrSend = mpz_get_str(NULL, 32, LocCurrFact);
        MPI_Request myRequest;
        MPI_Isend(FactStrSend, strlen(FactStrSend) + 1, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD, &myRequest);
    }

    

    // By now all processes have value of their maximum factorial stored in LocCurrFact
    // Convert all integer sums to floating point ones and perform division by largest factorial
    // to get true sum of THIS process
    // Precision chosen to be 64 + [ln(10)/ln(2) * N] bits
    mpf_set_default_prec(64 + ceil(3.33 * N));

    mpf_t LocSum_float;
    mpf_init(LocSum_float);
    mpf_set_z(LocSum_float, LocSum);

    mpf_t RankMaxFact_float;
    mpf_init(RankMaxFact_float);
    mpf_set_z(RankMaxFact_float, LocCurrFact);


    mpf_div(LocSum_float, LocSum_float, RankMaxFact_float);
    
    mpz_clears(LocCurrFact, LocSum, NULL);    
    mpf_clear(RankMaxFact_float);


    if(rank) {
        MPI_Send(&LocSum_float->_mp_prec, 1, MPI_INT,  0, 0, MPI_COMM_WORLD);
        MPI_Send(&LocSum_float->_mp_size, 1, MPI_INT,  0, 0, MPI_COMM_WORLD);
        MPI_Send(&LocSum_float->_mp_exp,  1, MPI_LONG, 0, 0, MPI_COMM_WORLD);


        int tmpLimbsSize = sizeof(LocSum_float->_mp_d[0]) * LocSum_float->_mp_size;

        MPI_Send(&tmpLimbsSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

        MPI_Send((char*)(LocSum_float->_mp_d), tmpLimbsSize, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

        mpf_clear(LocSum_float);
    }
    else {
        mpf_add_ui(LocSum_float, LocSum_float, 1);

        mpf_t Sum_i;        
        for(int i = 1; i < commsize; i++) {

            MPI_Recv(&Sum_i->_mp_prec, 1, MPI_INT,  i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&Sum_i->_mp_size, 1, MPI_INT,  i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&Sum_i->_mp_exp,  1, MPI_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


            int tmpSize;
            
            MPI_Recv(&tmpSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            mp_limb_t* tmpLimbs = (mp_limb_t*)calloc(tmpSize, sizeof(char));
    
            MPI_Recv((char*)(tmpLimbs), tmpSize, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
            Sum_i->_mp_d = tmpLimbs;       

            mpf_add(LocSum_float, LocSum_float, Sum_i);
            free(tmpLimbs);
        }
        
        char* formatStr = (char*)calloc(14 + strlen(argv[1]), sizeof(char));

        snprintf(formatStr, 13 + strlen(argv[1]), "%%.%dFf\b \b\n", N + 1);
        gmp_printf(formatStr, LocSum_float);

        free(formatStr);
        mpf_clear(LocSum_float);
    }

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}