#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from sys import argv

def main():
    if len(argv) < 3:
        print(f"Usage: {argv[0]} <orig_filename> <fft_filename>")
        return

    orig_signal = []
    fft_image = []

    with open(argv[1], 'r') as file:
        N = int(file.readline())
        
        for i in range(N):
            real = float(file.readline())
            imag = float(file.readline())
            orig_signal.append(complex(real, imag))

        orig_signal = np.array(orig_signal)

    with open(argv[2], 'r') as file:
        N = int(file.readline())
        
        for i in range(N):
            real = float(file.readline())
            imag = float(file.readline())
            fft_image.append(complex(real, imag))

        fft_image = np.array(fft_image)


    true_fft_image = np.fft.fft(orig_signal)

    if np.allclose(true_fft_image, fft_image, atol=1e-2):
        print('FFT is correct')
    else:
        print('FFT is not correct')
        print(np.max(np.abs(true_fft_image - fft_image)))



if __name__ == "__main__":
    main()
