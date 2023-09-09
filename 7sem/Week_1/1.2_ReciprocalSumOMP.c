#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define NUM_THREADS 4


int main(int argc, char* argv[])
{
    if(argc < 2) {
        printf("Usage: %s [N]\n", argv[0]);
        return 0;
    }

    double recip_sum = 0.0;

    // Passing number which we should sum up to by argv[1]
    const long long N = atoll(argv[1]);
    omp_set_dynamic(0);
    omp_set_num_threads(NUM_THREADS);

    #pragma omp parallel for reduction(+:recip_sum)
    for(long long i = 1; i < N + 1; i++) {
        recip_sum += 1.0 / (double)i;
    }

    printf("\nTotal sum S = %.16lg\n", recip_sum);

    return 0;
}