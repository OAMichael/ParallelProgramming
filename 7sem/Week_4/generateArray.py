#!/usr/bin/python3

import numpy as np
import sys

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <out_filename> [array_size]")
        return

    arraySize = 0
    if len(sys.argv) >= 3:
        arraySize = int(sys.argv[2])
    else:
        arraySize = np.random.randint(low=1000000, high=10000000)
    
    array = np.random.randint(low=-arraySize, high=arraySize, size=arraySize)

    with open(sys.argv[1], 'w') as file:
        file.write(str(arraySize) + "\n")
        for elem in array:
            file.write(str(elem) + "\n")




if __name__ == "__main__":
    main()