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

    factorial_t MaxFact = CalculateMaxN(N);

    mpz_t currFact;
    mpz_init_set(currFact, MaxFact.n_fact);

    mpz_t S;
    mpz_init_set(S, MaxFact.n_fact);

    for(unsigned i = 1; i <= MaxFact.n; i++) {
        
        mpz_cdiv_q_ui(currFact, currFact, i);

        mpz_add(S, S, currFact);
    }


    mpf_set_default_prec(64 + 8 * N);

    mpf_t S_float;
    mpf_init(S_float);
    mpf_set_z(S_float, S);

    mpf_t LargeNum;
    mpf_init(LargeNum);
    mpf_set_z(LargeNum, MaxFact.n_fact);

    mpf_div(S_float, S_float, LargeNum);

    char* formatStr = (char*)malloc(6 + strlen(argv[1]) + 1);
    memset(formatStr, 0, 6 + strlen(argv[1]) + 1);

    snprintf(formatStr, 6 + strlen(argv[1]), "%%.%dFf\n", N);
    gmp_printf(formatStr, S_float);

    free(formatStr);





    /*
    char* tmp = mpz_get_str(NULL, 10, S);

    char *buf = (char*)malloc(strlen(tmp) + 2);
    memset(buf, 0, strlen(tmp) + 2);

    strncpy(buf, tmp, 1);
    int len = strlen(buf);
    strcpy(buf + len, ".");
    len++;
    strncpy(buf + len, tmp + 1, 10000);

    strcpy(tmp, buf);
     
    printf("%s\n", buf);
    free(buf);
    */    


    /*
    char trueE[10004];
    FILE *fptr;
    if ((fptr = fopen("e.dat", "r")) == NULL) {
        printf("Error! File cannot be opened.");
        exit(1);
    }

    fscanf(fptr, "%[^\n]", trueE);
    printf("Data from the file:\n%s\n", trueE);
    fclose(fptr);
    */
    

    mpz_clear(S);
    mpz_clear(currFact);
    mpf_clear(S_float);
    mpf_clear(LargeNum);
    
    mpz_clear(MaxFact.n_fact);
  


    // Finalizing MPI
    MPI_Finalize();

    return 0;
}