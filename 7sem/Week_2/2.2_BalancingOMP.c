#include <stdio.h>
#include <omp.h>


#define NUM_ITERS   65
#define NUM_THREADS 4


//#define SCHEDULE_STATIC
//#define SCHEDULE_DYNAMIC
//#define SCHEDULE_GUIDED
#define SCHEDULE_DEFAULT

#define SCHEDULE_CHUNK_SIZE 4



#ifdef SCHEDULE_STATIC
    #define PARALLEL_STATEMENT schedule(static, SCHEDULE_CHUNK_SIZE)
#elif defined(SCHEDULE_DYNAMIC)
    #define PARALLEL_STATEMENT schedule(dynamic, SCHEDULE_CHUNK_SIZE)
#elif defined(SCHEDULE_GUIDED)
    #define PARALLEL_STATEMENT schedule(guided, SCHEDULE_CHUNK_SIZE)
#elif defined(SCHEDULE_DEFAULT)
    #define PARALLEL_STATEMENT
#endif



int main(int argc, char* argv[])
{
    omp_set_dynamic(0);
    omp_set_num_threads(NUM_THREADS);

    printf("Thread Iteration\n");

    #pragma omp parallel for PARALLEL_STATEMENT
    for(int i = 0; i < NUM_ITERS; i++) {
        printf("%d %d\n", omp_get_thread_num(), i);
    }

    return 0;
}