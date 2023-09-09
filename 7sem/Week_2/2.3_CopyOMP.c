#include <stdio.h>
#include <omp.h>

#define NUM_THREADS 4


#define USE_COPYIN
#define USE_COPYPRIVATE


int var_0 = 0;
int var_1 = 0;

#pragma omp threadprivate(var_0)
#pragma omp threadprivate(var_1)


int main(int argc, char* argv[])
{
    omp_set_dynamic(0);
    omp_set_num_threads(NUM_THREADS);

    printf("copyin:\n");
    #pragma omp parallel copyin(var_0)
    {
        #pragma omp single
        {
            var_0 = 1;
        }
        #pragma omp barrier
        printf("Thread %d: var = %d\n", omp_get_thread_num(), var_0);
    }


    printf("\ncopyprivate:\n");
    #pragma omp parallel
    {
        #pragma omp single copyprivate(var_1)
        {
            var_1 = 1;
        }
        #pragma omp barrier
        printf("Thread %d: var = %d\n", omp_get_thread_num(), var_1);
    }

    return 0;
}