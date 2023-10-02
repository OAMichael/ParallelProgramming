#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <omp.h>

#include "4.1_Sorts.h"

//#define SWITCH_SORT

//#define MAKE_TEST

#ifdef MAKE_TEST

#define AVERAGING_ITERATIONS	10

#define NUM_THREADS_ITERATIONS 	8
#define REC_DEPTH_ITERATIONS 	8

int g_RecDepth = NUM_THREADS_ITERATIONS;
int g_NumThreads = REC_DEPTH_ITERATIONS;

double g_TimeData[REC_DEPTH_ITERATIONS][NUM_THREADS_ITERATIONS] = {};

#else	// ^^^ MAKE_TEST ^^^   vvv !MAKE_TEST vvv

#define NUM_THREADS 2
int g_RecDepth = 6;

#endif	// ^^^ !MAKE_TEST ^^^


static void merge(int* array, const int left, const int mid, const int right)
{
    const int subArrayOne = mid - left + 1;
    const int subArrayTwo = right - mid;
 
    int *leftArray  = (int*)malloc(subArrayOne * sizeof(int)),
        *rightArray = (int*)malloc(subArrayTwo * sizeof(int));
 
    for(int i = 0; i < subArrayOne; i++)
        leftArray[i] = array[left + i];
    
    for(int i = 0; i < subArrayTwo; i++)
        rightArray[i] = array[mid + 1 + i];
 
    int indexOfSubArrayOne = 0,
        indexOfSubArrayTwo = 0;

    int indexOfMergedArray = left;
 
    while(indexOfSubArrayOne < subArrayOne && indexOfSubArrayTwo < subArrayTwo) {
        if(leftArray[indexOfSubArrayOne] <= rightArray[indexOfSubArrayTwo]) {
            array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
            indexOfSubArrayOne++;
        }
        else {
            array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
            indexOfSubArrayTwo++;
        }
        indexOfMergedArray++;
    }
 
    while(indexOfSubArrayOne < subArrayOne) {
        array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
        indexOfSubArrayOne++;
        indexOfMergedArray++;
    }

    while(indexOfSubArrayTwo < subArrayTwo) {
        array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
        indexOfSubArrayTwo++;
        indexOfMergedArray++;
    }
    free(leftArray);
    free(rightArray);
}
 

static void mergeSortRecursive(int* array, const int begin, const int end, const int depth)
{

#ifdef SWITCH_SORT

	if(depth < g_RecDepth) {
	    if(begin >= end)
	        return;
	 
	    const int mid = begin + (end - begin) / 2;

	    #pragma omp parallel
	    {
	    	#pragma omp single
	    	{
				#pragma omp task
			    mergeSortRecursive(array, begin, mid, depth + 1);

				#pragma omp task
		    	mergeSortRecursive(array, mid + 1, end, depth + 1);
	    	}
	    }
	    
	    merge(array, begin, mid, end);
	}
	else {
		heapSort(array + begin, end - begin + 1);
	}

#else

    if(begin >= end)
        return;
 
    const int mid = begin + (end - begin) / 2;

    #pragma omp parallel if(depth < g_RecDepth)
    {
    	#pragma omp single
    	{
			#pragma omp task
		    mergeSortRecursive(array, begin, mid, depth + 1);

			#pragma omp task
	    	mergeSortRecursive(array, mid + 1, end, depth + 1);
    	}
    }
    
    merge(array, begin, mid, end);

#endif

}


static void mergeSort(int* array, const int begin, const int end) {
	mergeSortRecursive(array, begin, end, 0);
}



int main(int argc, char** argv) {

	if (argc < 3) {
		printf("Usage: %s <in_filename> <sorted_filename> [time_data_filename]\n", argv[0]);
		return -1;
	}

    char strSize[32];

	FILE* fileIn = fopen(argv[1], "r"); 
    fgets(strSize, 32, fileIn);
    
    int arraySize = atoi(strSize);
    int* arrayToSort = (int*)malloc(arraySize * sizeof(int));

    for(int i = 0; i < arraySize; i++) {
        fgets(strSize, 32, fileIn);   
        arrayToSort[i] = atoi(strSize);
    }
    fclose(fileIn);

    printf("Input array size = %d\n", arraySize);

    omp_set_dynamic(0);
    omp_set_nested(1);

#ifdef MAKE_TEST

    for (int i = 0; i < AVERAGING_ITERATIONS; i++) {
	    for (int currentDepth = 1; currentDepth <= REC_DEPTH_ITERATIONS; currentDepth++) {
	    	for (int currentThreads = 1; currentThreads <= NUM_THREADS_ITERATIONS; currentThreads++) {
	    		int* testArray = (int*)malloc(arraySize * sizeof(int));
	    		memcpy(testArray, arrayToSort, arraySize * sizeof(int));

	    		g_RecDepth = currentDepth;
	    		g_NumThreads = currentThreads;

			    omp_set_num_threads(g_NumThreads);

			    double start = omp_get_wtime();
			    // Actual sorting
			    mergeSort(testArray, 0, arraySize - 1);
			    double end = omp_get_wtime();

			    g_TimeData[currentDepth - 1][currentThreads - 1] += end - start;
			    free(testArray);
			}

		    int percentage = i * REC_DEPTH_ITERATIONS + currentDepth;
		    printf("Completed: %.2f%%\n", (100.0f * percentage) / (AVERAGING_ITERATIONS * REC_DEPTH_ITERATIONS));
		}
	}


	FILE* fileTimeData;
	if (argc >= 4)  {
		fileTimeData = fopen(argv[3], "w");
		fprintf(fileTimeData, "Depth Threads Time\n");
	}

	for (int i = 0; i < REC_DEPTH_ITERATIONS; i++) {
		for (int j = 0; j < NUM_THREADS_ITERATIONS; j++) {
			g_TimeData[i][j] /= AVERAGING_ITERATIONS;
    		printf("Elapsed time for (RecDepth=%d, NumThreads=%d): %lf seconds\n", i + 1, j + 1, g_TimeData[i][j]);

    		if (argc >= 4) {
    			fprintf(fileTimeData, "%d %d %lf\n", i + 1, j + 1, g_TimeData[i][j]);
    		}
		}
	}
	if (argc >= 4) {
		fclose(fileTimeData);
	}

#else

    omp_set_num_threads(NUM_THREADS);

    double start = omp_get_wtime();
    // Actual sorting
    mergeSort(arrayToSort, 0, arraySize - 1);
    double end = omp_get_wtime();

    printf("Elapsed time for (RecDepth=%d, NumThreads=%d): %lf\n", g_RecDepth, NUM_THREADS, end - start);


    FILE* fileOut = fopen(argv[2], "w");
    fprintf(fileOut, "%d\n", arraySize);

    for (int i = 0; i < arraySize; i++) {
    	fprintf(fileOut, "%d\n", arrayToSort[i]);
    }
    fclose(fileOut);

#endif

    free(arrayToSort);
	return 0;
}
