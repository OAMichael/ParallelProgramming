#ifndef SORTS_H
#define SORTS_H

#include <stdbool.h>


// ========================= Bubble Sort ========================= //

void bubbleSort(int* arr, int n)
{
	int tmp = 0;
    bool swapped = false;

    for (int i = 0; i < n - 1; i++) {
        swapped = false;
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
            	tmp = arr[j];
            	arr[j] = arr[j + 1];
            	arr[j + 1] = tmp;
                swapped = true;
            }
        }
 
        if (!swapped)
            break;
    }
}


// ========================== Shell Sort ========================= //

int shellSort(int* arr, int n)
{
    for (int gap = n / 2; gap > 0; gap /= 2)
    {
        for (int i = gap; i < n; i += 1)
        {
            int temp = arr[i];
  
            int j;            
            for (j = i; j >= gap && arr[j - gap] > temp; j -= gap)
                arr[j] = arr[j - gap];
              
            arr[j] = temp;
        }
    }
    return 0;
}


// ========================== Comb Sort ========================== //

int getNextGap(int gap)
{
    gap = (gap * 10) / 13;
    if (gap < 1)
        return 1;

    return gap;
}
 
void combSort(int* arr, int n)
{
    int gap = n;
 
    bool swapped = true;
 
    while (gap != 1 || swapped)
    {
        gap = getNextGap(gap);
 
        swapped = false;
 
        for (int i = 0; i < n - gap; i++)
        {
            if (arr[i] > arr[i + gap])
            {
            	int tmp = arr[i];
            	arr[i] = arr[i + gap];
            	arr[i + gap] = tmp;
                swapped = true;
            }
        }
    }
}


// ========================== Heap Sort ========================== //

void heapify(int* arr, int N, int i)
{
 
    int largest = i;
 
    int l = 2 * i + 1;
    int r = 2 * i + 2;
 
    if (l < N && arr[l] > arr[largest])
        largest = l;
 
    if (r < N && arr[r] > arr[largest])
        largest = r;
 
    if (largest != i) {
    	int tmp = arr[i];
    	arr[i] = arr[largest];
    	arr[largest] = tmp;
 
        heapify(arr, N, largest);
    }
}
 
void heapSort(int* arr, int N)
{
 
    for (int i = N / 2 - 1; i >= 0; i--)
        heapify(arr, N, i);
 
    for (int i = N - 1; i > 0; i--) {
 
        int tmp = arr[0];
        arr[0] = arr[i];
        arr[i] = tmp;
 
        heapify(arr, i, 0);
    }
}


#endif // SORTS_H
