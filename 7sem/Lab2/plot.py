#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt


executors                   = [1,        2,        3,        4,        5,        6,        7,        8]
times_of_exec_baselineOMP   = [0.471270, 0.243355, 0.165881, 0.139147, 0.142626, 0.124460, 0.110018, 0.0993587]
times_of_exec_baselineMPI   = [0.474932, 0.293600, 0.231112, 0.204298, 0.257333, 0.236882, 0.230818, 0.2138020]
times_of_exec_task1_I       = [0.618433, 0.441768, 0.253863, 0.358923, 0.349970, 0.462114, 0.509942, 0.5322480]
times_of_exec_task1_J       = [0.621302, 0.404258, 0.357670, 0.352227, 0.411353, 0.433244, 0.466363, 0.5125010]
times_of_exec_task2         = [0.495985, 0.272414, 0.179548, 0.146769, 0.169057, 0.147873, 0.140251, 0.1528330]
times_of_exec_task3         = [0.526371, 0.307816, 0.235224, 0.200615, 0.219965, 0.202148, 0.190927, 0.1779900]



def plotSpeedup(ax, times, label):
    speedup_of_executors      = []
    efficiency_of_executors   = []
    
    for i in range(0, len(times)):
        speedup_of_executors.append(times[0] / times[i])
        efficiency_of_executors.append(times[0] / times[i] / executors[i])

    ax[0].plot(executors[0:len(times)], speedup_of_executors,    label=label, linewidth=2.5)
    ax[1].plot(executors[0:len(times)], efficiency_of_executors, label=label, linewidth=2.5)


def main():

    fig, ax = plt.subplots(nrows=2, figsize=[12, 10]) 

    fig.suptitle(r"Speedup $S$ and Efficiency $E$ as functions of $n$", fontsize=20)

    plotSpeedup(ax, times_of_exec_baselineOMP, "Baseline OMP")
    plotSpeedup(ax, times_of_exec_baselineMPI, "Baseline MPI")
    plotSpeedup(ax, times_of_exec_task1_I, "Task 1 parallel i (MPI)")
    plotSpeedup(ax, times_of_exec_task1_J, "Task 1 parallel j (MPI)")
    plotSpeedup(ax, times_of_exec_task2, "Task 2 (OMP)")
    plotSpeedup(ax, times_of_exec_task3, "Task 3 (OMP)")
    
    ax[0].grid()
    ax[0].set_xlabel(r"Executors $n$", fontsize=18)
    ax[0].set_ylabel(r"Speedup $S$", fontsize=18)
    ax[0].legend(loc="best", fontsize=14)

    ax[1].grid()
    ax[1].set_xlabel(r"Executors $n$", fontsize=18)
    ax[1].set_ylabel(r"Efficiency $E$", fontsize=18)

    plt.show()



if __name__ == '__main__':
    main()