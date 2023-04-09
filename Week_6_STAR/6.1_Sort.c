#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define PRINT_ARRAY


static inline void printArray(int* arr, size_t size) {
    for(int i = 0; i < size; i++)
        printf("%d ", arr[i]);
    printf("\n");
    return;
}


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
 

static void mergeSort(int* array, const int begin, const int end)
{
    if(begin >= end)
        return;
 
    const int mid = begin + (end - begin) / 2;
    mergeSort(array, begin, mid);
    mergeSort(array, mid + 1, end);
    merge(array, begin, mid, end);
}




int main(int argc, char* argv[]) {

    // Initializing MPI
    MPI_Init(&argc, &argv);

    // Getting the rank and communicator size
    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(argc < 2) {
        if(!rank)
            printf("Usage: %s <filename>\n", argv[0]);

        // Finalizing MPI
        MPI_Finalize();
        return 1;
    }


    FILE* fp;
    int arraySize = 0;
    char strSize[16];

    if(!rank) {
        fp = fopen(argv[1], "r"); 
        fgets(strSize, 16, fp);
        
        arraySize = atoi(strSize);
    }

    if(commsize > 1) {
        MPI_Bcast(&arraySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    int* arrayToSort = (int*)malloc(arraySize * sizeof(int));

    if(!rank) {
        for(int i = 0; i < arraySize; i++) {
            fgets(strSize, 16, fp);   
            arrayToSort[i] = atoi(strSize);
        }
        fclose(fp);

#ifdef PRINT_ARRAY
        printf("Given array:\n");
        printArray(arrayToSort, arraySize);
#endif
    }

    if(commsize > 1) {
        const int diff  = arraySize / commsize;
        int start = diff * rank;
        int end   = diff * (rank + 1) - 1;

        if(arraySize % commsize) {
            if(rank < arraySize % commsize) {
                start += rank;
                end   += rank + 1;
            }
            else {
                start += (arraySize % commsize);
                end   += (arraySize % commsize);
            }
        }

        int localArraySize = end - start + 1;
        int* localArray = (int*)malloc(localArraySize * sizeof(int));

        int* recvcnts = (int*)malloc(commsize * sizeof(int));
        int* displs = (int*)malloc(commsize * sizeof(int));
        int thisProcessSendCount = localArraySize;
        MPI_Gather(&thisProcessSendCount, 1, MPI_INT, recvcnts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        displs[0] = 0;
        for(int i = 1; i < commsize; ++i) {
            displs[i] = displs[i - 1] + recvcnts[i - 1]; 
        }


        MPI_Scatterv(arrayToSort, recvcnts, displs, MPI_INT, localArray, localArraySize, MPI_INT, 0, MPI_COMM_WORLD);
        mergeSort(localArray, 0, end - start);
        MPI_Gatherv(localArray, localArraySize, MPI_INT, arrayToSort, recvcnts, displs, MPI_INT, 0, MPI_COMM_WORLD);

        
        if(!rank) {
            int curr_start = 0; 
            int curr_mid   = recvcnts[0] - 1;
            int curr_end   = recvcnts[0] + recvcnts[1] - 1;

            for(int i = 0; i < commsize - 1; i++) {
                merge(arrayToSort, curr_start, curr_mid, curr_end);
                curr_end += recvcnts[i + 1];
                curr_mid += recvcnts[i + 1];
            }
        }


        free(recvcnts);
        free(displs);
        free(localArray);
    }
    else {
        mergeSort(arrayToSort, 0, arraySize - 1);
    }

#ifdef PRINT_ARRAY
    if(!rank) {        
        printf("Sorted array:\n");
        printArray(arrayToSort, arraySize);
    }
#endif

    free(arrayToSort);

    // Finalizing MPI
    MPI_Finalize();

    return 0;
}