#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


speedup_effiency_plot = False

processes    = [1, 		2, 		3, 		4, 		5, 		6, 		7, 		8]
times        = []

arrayN       = [1000000, 10000000, 100000000]

# 1000000
times.append([0.243088, 0.141104, 0.106997, 0.0946558, 0.0986454, 0.103553, 0.102976, 0.104481])

# 10000000
times.append([2.78789, 1.68819, 1.24211, 1.12385, 1.21597, 1.46908, 1.45649, 1.41982])

# 100000000
times.append([34.3703, 22.6197, 18.5721, 18.7826, 18.1110, 18.072, 18.0169, 18.8009])

speedup      = []
efficiency   = []

def main():

	fig, ax = plt.subplots(nrows=2, figsize=[12, 10]) 

	fig.suptitle(r"Speedup $S$ and Efficiency $E$ as functions of $p$", fontsize=20)


	for i in range(0, len(times)):
		new_speedup    = []
		new_efficiency = []
		for j in range(0, len(processes)):
			new_speedup.append(times[i][0] / times[i][j])
			new_efficiency.append(times[i][0] / times[i][j] / processes[j])

		speedup.append(new_speedup)
		efficiency.append(new_efficiency)

		ax[0].plot(processes, new_speedup,    label=f"n = {arrayN[i]}")
		ax[1].plot(processes, new_efficiency, label=f"n = {arrayN[i]}")

	ax[0].grid()
	ax[0].set_xlabel(r"Processes $p$", fontsize=18)
	ax[0].set_ylabel(r"Speedup $S$", fontsize=18)
	ax[0].legend(loc="upper left", fontsize=14)


	ax[1].grid()
	ax[1].set_xlabel(r"Processes $p$", fontsize=18)
	ax[1].set_ylabel(r"Efficiency $E$", fontsize=18)

	plt.show()



if __name__ == '__main__':
    main()