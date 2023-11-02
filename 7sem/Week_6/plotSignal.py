#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from sys import argv

def main():
    if len(argv) < 2:
        print(f"Usage: {argv[0]} <in_filename>")
        return

    real = []
    imag = []

    with open(argv[1], 'r') as file:
        N = int(file.readline())

        for i in range(N):
            real.append(float(file.readline()))
            imag.append(float(file.readline()))

    fig, _ = plt.subplots(figsize=[16, 9])
    plt.stem(np.linspace(0, N, N), real, 'blue')
    plt.stem(np.linspace(0, N, N), imag, 'orange')
    plt.grid()
    plt.title('Sine wave')
    plt.show()
 


if __name__ == "__main__":
    main()
