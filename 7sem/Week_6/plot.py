#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt


time_fft_base = False
use_log_scale = False

fft_bases = np.array([ 4096,
                       8192,
                       16384,
                       32768,
                       65536,
                       131072,
                       262144,
                       524288,
                       1048576,
                       2097152 ])

times_of_fft_base_seq = np.array([ 0.007314,
                                   0.010947,
                                   0.023747,
                                   0.049073,
                                   0.102170,
                                   0.213532,
                                   0.445707,
                                   0.941815,
                                   1.954418,
                                   4.011739 ])

times_of_fft_base_seq_vec = np.array([ 0.003398,
                                       0.008541,
                                       0.015998,
                                       0.033061,
                                       0.069762,
                                       0.146413,
                                       0.309494,
                                       0.703970,
                                       1.370407,
                                       2.865753 ])
                            
times_of_fft_base_par = np.array([ 0.005174,
                                   0.010899,
                                   0.018627,
                                   0.029062,
                                   0.060978,
                                   0.120441,
                                   0.238550,
                                   0.500154,
                                   1.028081,
                                   2.125929 ])

times_of_fft_base_par_vec = np.array([ 0.003732,
                                       0.007146,
                                       0.011949,
                                       0.020753,
                                       0.042439,
                                       0.087030,
                                       0.178966,
                                       0.382204,
                                       0.787187,
                                       1.629634 ])



speedup_effiency_plot = True

threads                 = [1,         2,        3,        4,        5,        6,        7,        8]
times_of_threads        = [1.925180,  1.024923, 1.030226, 1.071849, 1.071883, 1.159222, 1.130189, 1.123179]
speedup_of_threads      = []
efficiency_of_threads   = []
fixed_fft_base          = 1048576

def main():

    if time_fft_base:

        fig = plt.figure(figsize=[16, 8])
        ax = fig.add_subplot(111)

        fig.suptitle(r"Elapsed time for FFT (2 threads for parallel versions)", fontsize=24)

        ax.grid(which='major', color="black", linewidth=1.5)
        ax.grid(which='minor', color="grey",  linewidth=0.2)
        
        plt.ticklabel_format(style='plain') 
        
        if use_log_scale:
            ax.set_xscale("log")
            ax.set_yscale("log")

        ax.plot(fft_bases, times_of_fft_base_seq, linewidth=2.5, label='Sequential')
        ax.plot(fft_bases, times_of_fft_base_seq_vec, linewidth=2.5, label='Sequential, vector instructions')
        ax.plot(fft_bases, times_of_fft_base_par, linewidth=2.5, label='Parallel')
        ax.plot(fft_bases, times_of_fft_base_par_vec, linewidth=2.5, label='Parallel, vector instructions')
    
        ax.set_xlabel(r"$n$", fontsize=20)
        ax.set_ylabel(r"$t$, sec", fontsize=20)
        ax.autoscale()

        plt.legend(loc='best', fontsize=18)
        plt.show()

    if speedup_effiency_plot:
        fig, ax = plt.subplots(nrows=2, figsize=[12, 10]) 

        fig.suptitle(r"Speedup $S$ and Efficiency $E$ as functions of $t$", fontsize=20)


        for i in range(0, len(threads)):
            speedup_of_threads.append(times_of_threads[0] / times_of_threads[i])
            efficiency_of_threads.append(times_of_threads[0] / times_of_threads[i] / threads[i])

        ax[0].plot(threads, speedup_of_threads,    label=f"n = {fixed_fft_base}")
        ax[1].plot(threads, efficiency_of_threads, label=f"n = {fixed_fft_base}")

        ax[0].grid()
        ax[0].set_xlabel(r"Threads $t$", fontsize=18)
        ax[0].set_ylabel(r"Speedup $S$", fontsize=18)
        ax[0].legend(loc="best", fontsize=14)

        ax[1].grid()
        ax[1].set_xlabel(r"Threads $t$", fontsize=18)
        ax[1].set_ylabel(r"Efficiency $E$", fontsize=18)

        plt.show()



if __name__ == '__main__':
    main()