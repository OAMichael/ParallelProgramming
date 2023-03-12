#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt, ticker as mticker
import pandas as pd



def main():
	df = pd.read_csv("../build/Data.txt", sep=" ")

	size = df.iloc[:, 0]
	time = df.iloc[:, 1]

	time = time * (10 ** 6)
	fig, ax = plt.subplots(figsize=[12, 8])

	fig.suptitle("Send time is a variable of send size", fontsize=20)
	ax.plot(size, time)
	ax.set_xlabel("Send size in bytes", fontsize=16)
	ax.set_ylabel(r"Send time, $\mu s$", fontsize=16)

	ax.set_xscale('log')
	ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
	ax.grid()

	plt.show()



if __name__ == '__main__':
    main()