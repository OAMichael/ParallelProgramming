#!/usr/bin/python3

import numpy as np
import sys


def main():

    if len(sys.argv) < 4:
        print(f"Usage: {argv[0]} <A_filename> <B_filename> <C_filename>")
        return

    A = []
    B = []
    C = []

    with open(sys.argv[1], 'r') as f:
        A = [[float(num) for num in line.split()] for line in f]

    with open(sys.argv[2], 'r') as f:
        B = [[float(num) for num in line.split()] for line in f]

    with open(sys.argv[3], 'r') as f:
        C = [[float(num) for num in line.split()] for line in f]

    A = np.array(A)
    B = np.array(B)
    C = np.array(C)

    print("A:")
    print(A)
    print("B:")
    print(B)

    print("\nC:")
    print(C)


    trueResult = np.matmul(A, B)
    if np.allclose(C, trueResult, 1e-3):
        print("Result is correct")
    else:
        print("Result is incorrect")
        print(C - trueResult)




if __name__ == "__main__":
    main()