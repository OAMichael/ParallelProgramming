#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>
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
    mpz_set_ui(n, (unsigned int)(ull >> 32)); /* n = (unsigned int)(ull >> 32) */
    mpz_mul_2exp(n, n, 32);                   /* n <<= 32 */
    mpz_add_ui(n, n, (unsigned int)ull);      /* n += (unsigned int)ull */
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

    //int_fast64_t a = end - 1;
    //int_fast64_t LocSumUint = end;

    mpz_t a_mpz;
    mpz_init_set_ui(a_mpz, end - 1);

    mpz_t S_mpz;
    mpz_init_set_ui(S_mpz, end);

    // Main algorithm
    mpz_t LocCurrFact;
    mpz_init_set_ui(LocCurrFact, 1);

    mpz_t LocSum;
    //mpz_init_set_ui(LocSum, end);
    //mpz_init(LocSum);
    mpz_init_set_ui(LocSum, 0);

    

    
    //clock_t start_time = clock();
    

    for(int i = end - 2; i > start; i--) {
        if(/*a > UINT64_MAX / (int_fast64_t)i */ i % 8192 == 0) {
            //mpz_addmul_ui(LocSum, LocCurrFact, LocSumUint);
            //mpz_mul_ui(LocCurrFact, LocCurrFact, a);

            //a = (int_fast64_t)i;
            //LocSumUint = a;

            mpz_addmul(LocSum, LocCurrFact, S_mpz);
            mpz_mul(LocCurrFact, LocCurrFact, a_mpz);

            mpz_set_ull(a_mpz, i);
            mpz_set(S_mpz, a_mpz);
        }
        else {
            mpz_mul_ui(a_mpz, a_mpz, i);
            mpz_add(S_mpz, S_mpz, a_mpz);
            //a *= (int_fast64_t)i;
            //LocSumUint += a;
        }
    }

    //mpz_addmul_ui(LocSum, LocCurrFact, LocSumUint);
    //mpz_mul_ui(LocCurrFact, LocCurrFact, a);
    mpz_addmul(LocSum, LocCurrFact, S_mpz);
    mpz_mul(LocCurrFact, LocCurrFact, a_mpz);
    
    //clock_t end_time = clock();
    //printf("Process %d: main loop time elapsed: %lf\n", rank, ((double)(end_time - start_time))/CLOCKS_PER_SEC);

    //MPI_Finalize();
    //return 0;
    


    /*
    for(int i = end - 2; i > start; i--) {
        if(a > UINT64_MAX / (uint64_t)i) {
            mpz_mul_ui(LocCurrFact, LocCurrFact, a);
            mpz_addmul_ui(LocSum, LocCurrFact, i);
            a = (uint64_t)i;
        }
        else {
            a *= (uint64_t)i;
            mpz_addmul_ui(LocSum, LocCurrFact, a);
        }
        mpz_mul_ui(LocCurrFact, LocCurrFact, i);
    }
    mpz_mul_ui(LocCurrFact, LocCurrFact, a);
    */
    
    /*
    for(int i = end - 2; i > start; i--) {
        if(!isOverflow) {
            if(a > UINT64_MAX / (uint64_t)i) {
                isOverflow = true;
                mpz_set_ui(LocCurrFact, a);
                //mpz_set_ui(LocSum, LocSumUint);

                mpz_mul_ui(LocCurrFact, LocCurrFact, i);
                mpz_add(LocSum, LocSum, LocCurrFact);


                //mpz_addmul_ui(LocSum, LocCurrFact, i);
            }   
            else {
                a *= (uint64_t)i;
                //LocSumUint += a;
                mpz_add_ui(LocSum, LocSum, a);
            }
        }
        else {
            mpz_mul_ui(LocCurrFact, LocCurrFact, i);
            mpz_add(LocSum, LocSum, LocCurrFact);
        }
    }
    */

    //if(!isOverflow)
        //mpz_set_ui(LocCurrFact, a);





    // For all we calculate LocCurrFact = start * (start + 1) * ... * (end - 2) * (end - 1)
    mpz_mul_ui(LocCurrFact, LocCurrFact, start);


    
    //clock_t start_time = clock();
    

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

    
    //clock_t end_time = clock();
    //printf("Process %d: main loop time elapsed: %lf\n", rank, ((double)(end_time - start_time))/CLOCKS_PER_SEC);

    //MPI_Finalize();
    //return 0;
    
    
    //clock_t start_time = clock();
    


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
    //clock_t end_time = clock();
    //printf("Process %d: main loop time elapsed: %lf\n", rank, ((double)(end_time - start_time))/CLOCKS_PER_SEC);

    //MPI_Finalize();
    //return 0;
    mpz_clear(LocCurrFact);
    mpz_clear(LocSum);
    mpf_clear(RankMaxFact_float);
    


    //clock_t start_time = clock();
    /*
    // Send sums from all processes to 0's process
    if(rank) {

        char* buf = (char*)calloc(N + 5, sizeof(char));

        char* formatStr = (char*)calloc(5 + strlen(argv[1]) + 1, sizeof(char));
        snprintf(formatStr, 5 + strlen(argv[1]) + 1, "%%.%dFf", N + 2);

        gmp_snprintf(buf, N + 5, formatStr, LocSum_float);
        free(formatStr);

        MPI_Send(buf, strlen(buf) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        free(buf);
    }
    if(!rank) {
        // Add 1 to first process because we need to take into account 0! = 1
        mpf_add_ui(LocSum_float, LocSum_float, 1);

        // Now accumulate all sums from all processes
        // Still valid if we have only 1 process
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

        
        // Make nice looking format string. strlen(argv[1]) + 1 because N can be 999, so N + 1 = 1000
        // Allocate [14 + strlen(argv[1])] bytes to contain format string: 
        // %%.      -- 3 bytes
        // %d       -- itoa(N + 1) <= strlen(argv[1]) + 1
        // Ff       -- 2 bytes
        // \b \b\n  -- 2 + 1 + 2 + 2 = 7 bytes
        // Total: <= 13 + strlen(argv[1])
        // +1 for '\0'
        char* formatStr = (char*)calloc(14 + strlen(argv[1]), sizeof(char));

        snprintf(formatStr, 13 + strlen(argv[1]), "%%.%dFf\b \b\n", N + 1);
        gmp_printf(formatStr, LocSum_float);

        free(formatStr);
        mpf_clear(LocSum_float);
    }
    */

    if(rank) {
        mp_exp_t SendExp;
        char* SendBuf = mpf_get_str(NULL, &SendExp, 32, ceil(N * 0.7), LocSum_float);
        SendExp *= -1;
        MPI_Send(&SendExp, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(SendBuf, strlen(SendBuf) + 1 - SendExp, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    else {
        mpf_add_ui(LocSum_float, LocSum_float, 1);

        mp_exp_t RecvExp;

        mpf_t Sum_i;
        mpf_init(Sum_i);

        int recvLength = ceil(N * 0.7) + 3;
        char* strSum_i = (char*)calloc(recvLength + 1, sizeof(char));
        strSum_i[0] = '0';
        strSum_i[1] = '.';

        for(int i = 1; i < commsize; i++) {
            MPI_Recv(&RecvExp, 1, MPI_LONG, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int j = 2; j < 2 + RecvExp; j++)
                strSum_i[j] = '0';

            MPI_Recv(strSum_i + 2 + RecvExp, recvLength - RecvExp, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            mpf_set_str(Sum_i, strSum_i, 32);
            mpf_add(LocSum_float, LocSum_float, Sum_i);
        }
        free(strSum_i);

        mpf_clear(Sum_i);

        char* formatStr = (char*)calloc(14 + strlen(argv[1]), sizeof(char));

        snprintf(formatStr, 13 + strlen(argv[1]), "%%.%dFf\b \b\n", N + 1);
        gmp_printf(formatStr, LocSum_float);

        free(formatStr);
        mpf_clear(LocSum_float);
    }

    //clock_t end_time = clock();
    //printf("Process %d: main loop time elapsed: %lf\n", rank, ((double)(end_time - start_time))/CLOCKS_PER_SEC);


    // Finalizing MPI
    MPI_Finalize();

    return 0;
}