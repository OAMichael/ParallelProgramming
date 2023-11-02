#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from sys import argv

def main():
    if len(argv) < 3:
        print(f"Usage: {argv[0]} <N> <out_filename>")
        return

    N = int(argv[1])
    time = np.linspace(0, N, N)
    signal = np.sin(time * 2 * np.pi / 512)

    with open(argv[2], 'w') as file:
        signal_size = len(signal)

        file.write(str(signal_size) + "\n")
        for i in range(signal_size):
            file.write(str(signal[i].real) + "\n")
            file.write(str(signal[i].imag) + "\n")


if __name__ == "__main__":
    main()
