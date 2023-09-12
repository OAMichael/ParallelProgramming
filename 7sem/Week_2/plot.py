#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import sys


def main():
	filename = sys.argv[1]
	chunksize = int(sys.argv[2])
	title = sys.argv[3]

	df = pd.read_csv(filename, sep=" ")
	df.sort_values(by=['Iteration'], inplace=True)

	thread    = np.array(df.iloc[:, 0].tolist())
	iteration = np.array(df.iloc[:, 1].tolist())

	fig = plt.figure(figsize=[12, 8])
	fig.suptitle(f"{title} (chunk = {chunksize})")
	ax = plt.axes()

	im = ax.scatter(iteration, thread)
	ax.set_xlabel(r"Iteration", fontsize=18)
	ax.set_ylabel(r"Thread number", fontsize=18)
	ax.grid()
	ax.set_xlim([-1, 65])
	ax.set_ylim([-1, 5])

	plt.show()


if __name__ == '__main__':
    main()