#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


speedup_effiency_plot = False

processes    = [1, 	2, 	3, 	4, 	5, 	6, 	7,  8]
times        = []

end_points  = [r"$1-10^{-9}$", r"$1-10^{-9}$", r"$1-10^{-13}$"]
epsilons  	= [r"$10^{-3}$", r"$10^{-6}$", r"$10^{-3}$"]

# 0.99999999, 10^{-3}
times.append([4.096, 2.431, 1.648, 1.322, 1.206, 1.165, 1.156, 1.100])

# 0.99999999, 10^{-6}
times.append([5.094, 2.962, 2.005, 1.645, 1.444, 1.354, 1.269, 1.255])

# 0.999999999999, 10^{-3}
times.append([5.634, 3.227, 2.205, 1.976, 1.887, 1.768, 1.674, 1.399])

speedup      = []
efficiency   = []

def main():

	fig, ax = plt.subplots(nrows=2, figsize=[12, 10]) 

	fig.suptitle(r"Speedup $S$ and Efficiency $E$ as functions of number of threads $t$", fontsize=20)


	for i in range(0, len(times)):
		new_speedup    = []
		new_efficiency = []
		for j in range(0, len(processes)):
			new_speedup.append(times[i][0] / times[i][j])
			new_efficiency.append(times[i][0] / times[i][j] / processes[j])

		speedup.append(new_speedup)
		efficiency.append(new_efficiency)

		ax[0].plot(processes, new_speedup,    label=f"end = {end_points[i]}, " + r"$\varepsilon = $" + f"{epsilons[i]}")
		ax[1].plot(processes, new_efficiency, label=f"end = {end_points[i]}, " + r"$\varepsilon = $" + f"{epsilons[i]}")



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