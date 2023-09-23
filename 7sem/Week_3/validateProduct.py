#!/usr/bin/python3

import numpy as np


def main():

    A = []
    B = []
    C = []

    with open('../../build/7sem/3.1_data/A.dat', 'r') as f:
        A = [[float(num) for num in line.split()] for line in f]

    with open('../../build/7sem/3.1_data/B.dat', 'r') as f:
        B = [[float(num) for num in line.split()] for line in f]

    with open('../../build/7sem/3.1_data/C.dat', 'r') as f:
        C = [[float(num) for num in line.split()] for line in f]

    A = np.array(A)
    B = np.array(B)
    C = np.array(C)

    trueResult = np.matmul(A, B)
    if np.allclose(C, trueResult):
        print("Result is correct")
    else:
        print("Result is incorrect")
        print(C - trueResult)




if __name__ == "__main__":
    main()