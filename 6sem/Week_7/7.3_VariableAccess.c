#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <semaphore.h>

#define NUM_THREADS 4


typedef struct 
{
    sem_t* p_semArray;
    int* p_var;
    int tid;

} thread_info_t;


void* increment_var(void* arg) {

    thread_info_t* local_info = (thread_info_t*)arg;

    sem_wait(&(local_info->p_semArray[local_info->tid]));

    printf("Thread %d before increment: var = %d\n", local_info->tid, *local_info->p_var);
    *(local_info->p_var) += 1;
    printf("Thread %d after increment:  var = %d\n\n",  local_info->tid, *local_info->p_var);

    if(local_info->tid < NUM_THREADS - 1)
        sem_post(&(local_info->p_semArray[local_info->tid + 1]));

    return NULL;
}



int main(int argc, char* argv[])
{
    pthread_t threads[NUM_THREADS];
    thread_info_t infos[NUM_THREADS];
    sem_t semArray[NUM_THREADS];

    int var_to_increment = 0;

    for(int i = 0; i < NUM_THREADS; i++)
    {
        sem_init(&semArray[i], 0, i > 0 ? 0 : 1);

        infos[i].p_var = &var_to_increment;
        infos[i].tid = i;
        infos[i].p_semArray = semArray;
    }

    int rc = 0;
    for(int i = 0; i < NUM_THREADS; i++) {
        if((rc = pthread_create(&threads[i], NULL, increment_var, &infos[i]))) {
            perror("pthread_create");
            return -1;
        }
    }

    for(int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);

        sem_destroy(&semArray[i]);
    }



    return 0;
}