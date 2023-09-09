#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_THREADS 4


typedef struct 
{
    pthread_mutex_t* p_mutex;

    long long start, end;
    double* p_total_sum;

} thread_info_t;


// Function to sum reciprocals
void* recip_sum(void* arg) {

    thread_info_t* local_info = (thread_info_t*)arg;

    long long local_start = local_info->start;
    long long local_end   = local_info->end;

    if(local_start >= local_end)
        return NULL;

    double local_sum = 0.0;

    for(long long i = local_start; i < local_end; i++) {
        local_sum += 1.0 / (double)i;
    }

    pthread_mutex_lock(local_info->p_mutex);
    *(local_info->p_total_sum) += local_sum;
    pthread_mutex_unlock(local_info->p_mutex);

    return NULL;
}



int main(int argc, char* argv[])
{
    if(argc < 2) {
        printf("Usage: %s [N]\n", argv[0]);
        return 0;
    }

    // Initialize mutex
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_t tids[NUM_THREADS];

    double total_sum = 0.0;

    thread_info_t infos[NUM_THREADS];

    // Passing number which we should sum up to by argv[1]
    const long long N = atoll(argv[1]);

    // Distributing starts and ends of summing among processes
    const long long diff  = N / NUM_THREADS;

    for(int i = 0; i < NUM_THREADS; i++)
    {
        infos[i].p_mutex = &mutex;
        infos[i].p_total_sum = &total_sum;

        infos[i].start    = 1 + i * diff;
        infos[i].end      = 1 + (i + 1) * diff;

        if(N % NUM_THREADS) {
            if(i < N % NUM_THREADS) {
                infos[i].start += i;
                infos[i].end   += i + 1;
            }
            else {
                infos[i].start += N % NUM_THREADS;
                infos[i].end   += N % NUM_THREADS;
            }
        }
    }

    int rc = 0;
    for(int i = 0; i < NUM_THREADS; i++) {
        if((rc = pthread_create(&tids[i], NULL, recip_sum, &infos[i]))) {
            perror("pthread_create");
            return -1;
        }
    }

    for(int i = 0; i < NUM_THREADS; i++)
        pthread_join(tids[i], NULL);

    printf("\nTotal sum S = %.16lg\n", total_sum);

    return 0;
}