#include <stdio.h>
#include <omp.h>

#define NEST_0_NUM_THREADS 8
#define NEST_1_NUM_THREADS 2
#define NEST_2_NUM_THREADS 4


int main(int argc, char* argv[])
{
    omp_set_dynamic(0);
    omp_set_nested(1);

    printf("<current nest thread number>/<current nest number of threads>/[<previous nest number of threads>/...]\n\n");
    int sync_var = 0;

    printf("Nest level: 0\n");
    #pragma omp parallel num_threads(NEST_0_NUM_THREADS) shared(sync_var)
    {
        printf("%d/%d\n", omp_get_thread_num(), omp_get_num_threads());
        #pragma omp barrier

        if(omp_get_thread_num() == 0) {
            printf("Nest level: 1\n");
        }
        int nest0_num_threads = omp_get_num_threads();
        #pragma omp barrier

        #pragma omp parallel num_threads(NEST_1_NUM_THREADS) shared(nest0_num_threads)
        {    
            printf("%d/%d/%d\n", omp_get_thread_num(), omp_get_num_threads(), nest0_num_threads);
            #pragma omp barrier

            #pragma omp atomic
            sync_var++;
            #pragma omp flush(sync_var)
            while(sync_var < NEST_0_NUM_THREADS * NEST_1_NUM_THREADS) {
                #pragma omp flush(sync_var)
            }
            
            if(omp_get_thread_num() == 0  && omp_get_ancestor_thread_num(omp_get_level() - 1) == 0) {
                printf("Nest level: 2\n");
            }
            int nest1_num_threads = omp_get_num_threads();
            #pragma omp barrier

            #pragma omp parallel num_threads(NEST_2_NUM_THREADS) shared(nest1_num_threads)
            {
                printf("%d/%d/%d/%d\n", omp_get_thread_num(), omp_get_num_threads(), nest1_num_threads, nest0_num_threads);
            }
        }
    }
    return 0;
}