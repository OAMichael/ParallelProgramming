#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_THREADS 4


typedef struct 
{
    pthread_mutex_t* p_mutex;
    int* p_var;
    int tid;

} thread_info_t;


// Function to sum reciprocals
void* increment_var(void* arg) {

    thread_info_t* local_info = (thread_info_t*)arg;

    pthread_mutex_lock(local_info->p_mutex);

    printf("Thread %d before increment: var = %d\n", local_info->tid, *local_info->p_var);
    *(local_info->p_var) += 1;
    printf("Thread %d after increment:  var = %d\n\n",  local_info->tid, *local_info->p_var);

    pthread_mutex_unlock(local_info->p_mutex);

    return NULL;
}



int main(int argc, char* argv[])
{
    // Initialize mutex
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_t tids[NUM_THREADS];

    thread_info_t infos[NUM_THREADS];

    int var_to_increment = 0;

    for(int i = 0; i < NUM_THREADS; i++)
    {
        infos[i].p_mutex = &mutex;
        infos[i].p_var = &var_to_increment;
        infos[i].tid = i;
    }

    int rc = 0;
    for(int i = 0; i < NUM_THREADS; i++) {
        if((rc = pthread_create(&tids[i], NULL, increment_var, &infos[i]))) {
            perror("pthread_create");
            return -1;
        }
    }

    for(int i = 0; i < NUM_THREADS; i++)
        pthread_join(tids[i], NULL);

    return 0;
}