#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt


time_of_matrix_size = False
use_log_scale = False

matrix_sizes = np.array([ 1,
                          2,
                          4,
                          8,
                          16,
                          32,
                          64,
                          128,
                          256,
                          512,
                          1024,
                          2048,
                          4096,
                          8192 ])

times_of_size = np.array([ 0.000212,
                           0.000352,
                           0.000212,
                           0.000205,
                           0.000211,
                           0.000613,
                           0.002640,
                           0.009537,
                           0.026146,
                           0.106804,
                           0.727148,
                           5.650605,
                           47.109121,
                           371.260483 ])



speedup_effiency_plot = True

threads                 = [1,                    2,         3,         4,         5,         6,         7,         8]
times_of_threads        = [156.155347,  100.504811, 73.380066, 50.906166, 56.581669, 57.766249, 48.279266, 49.690378]
speedup_of_threads      = []
efficiency_of_threads   = []
fixed_matrix_size       = 4096

def main():

    if time_of_matrix_size:

        fig = plt.figure(figsize=[12, 8])
        ax = fig.add_subplot(111)

        fig.suptitle(r"Elapsed time of matrix multiplication, 4 threads", fontsize=24)

        ax.grid(which='major', color="black", linewidth=1.5)
        ax.grid(which='minor', color="grey", linewidth=0.2)
        ax.plot(matrix_sizes, times_of_size, linewidth=2.5)

        if use_log_scale:
            ax.set_xscale("log")
            ax.set_yscale("log")
    
        ax.set_xlabel(r"$n$", fontsize=20)
        ax.set_ylabel(r"$t$, sec", fontsize=20)
        ax.autoscale()

        plt.show()

    if speedup_effiency_plot:
        fig, ax = plt.subplots(nrows=2, figsize=[12, 10]) 

        fig.suptitle(r"Speedup $S$ and Efficiency $E$ as functions of $t$", fontsize=20)


        for i in range(0, len(threads)):
            speedup_of_threads.append(times_of_threads[0] / times_of_threads[i])
            efficiency_of_threads.append(times_of_threads[0] / times_of_threads[i] / threads[i])

        ax[0].plot(threads, speedup_of_threads,    label=f"n = {fixed_matrix_size}")
        ax[1].plot(threads, efficiency_of_threads, label=f"n = {fixed_matrix_size}")

        ax[0].grid()
        ax[0].set_xlabel(r"Threads $t$", fontsize=18)
        ax[0].set_ylabel(r"Speedup $S$", fontsize=18)
        ax[0].legend(loc="upper left", fontsize=14)

        ax[1].grid()
        ax[1].set_xlabel(r"Threads $t$", fontsize=18)
        ax[1].set_ylabel(r"Efficiency $E$", fontsize=18)

        plt.show()



if __name__ == '__main__':
    main()