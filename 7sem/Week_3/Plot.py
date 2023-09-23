#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


def main():

    df = pd.read_csv("../../build/7sem/3.1_data/4threads.dat", sep=" ")

    n  = np.array(df.iloc[:, 0].tolist())
    m1 = np.array(df.iloc[:, 1].tolist())
    m2 = np.array(df.iloc[:, 2].tolist())
    m3 = np.array(df.iloc[:, 3].tolist())
    m4 = np.array(df.iloc[:, 4].tolist())
    m5 = np.array(df.iloc[:, 5].tolist())
    m6 = np.array(df.iloc[:, 6].tolist())
    m7 = np.array(df.iloc[:, 7].tolist())

    fig = plt.figure(figsize=[12, 8])

    plt.title("Matrix multiplication methods comparison", fontsize=28)

    plt.plot(n, m1, label='Naive', linewidth=2.5)
    plt.plot(n, m2, label='Transpose', linewidth=2.5)
    plt.plot(n, m3, label='Sum and precopy', linewidth=2.5)
    plt.plot(n, m4, label='Naive parallel', linewidth=2.5)
    plt.plot(n, m5, label='Transpose parallel', linewidth=2.5)
    plt.plot(n, m6, label='Sum and precopy parallel', linewidth=2.5)
    plt.plot(n, m7, label='Naive parallel collapse', linewidth=2.5)

    plt.legend(loc="upper left", fontsize=15)
    plt.grid()

    plt.xlabel(r"$n$", fontsize=24)
    plt.ylabel(r"$t, s$", fontsize=24)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    plt.show()



if __name__ == '__main__':
    main()