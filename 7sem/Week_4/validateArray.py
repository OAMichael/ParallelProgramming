#!/usr/bin/python3

import numpy as np
import sys

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <in_filename> <sorted_filename>")
        return

    initialArray = []
    sortedArray = []

    with open(sys.argv[1], 'r') as file:
        initialArraySize = int(file.readline())
        initialArray = np.zeros(initialArraySize, dtype=int)

        for i in range(0, initialArraySize):
            initialArray[i] = int(file.readline())


    with open(sys.argv[2], 'r') as file:
        sortedArraySize = int(file.readline())
        sortedArray = np.zeros(sortedArraySize, dtype=int)

        for i in range(0, sortedArraySize):
            sortedArray[i] = int(file.readline())


    if initialArraySize != sortedArraySize:
        print("Array sizes are not equal!")
        return


    trueSorted = np.sort(initialArray)

    if np.array_equal(trueSorted, sortedArray):
        print("Arrays are equal")
    else:
        print("Arrays are not equal")


if __name__ == "__main__":
    main()