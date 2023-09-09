#include <stdio.h>
#include <omp.h>

#define NUM_THREADS 4


int main(int argc, char** argv) {
    omp_set_dynamic(0);
    omp_set_num_threads(NUM_THREADS);

	#pragma omp parallel
	printf("Hello from thread %d/%d\n", omp_get_thread_num(), NUM_THREADS);

	return 0;
}