#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


speedup_effiency_plot = True

processes    = [1, 		2, 		3, 		4, 		5, 		6, 		7, 		8]
times        = []

mesh_ns      = [2000, 4000, 8000]

# 2000
times.append([6.1414, 2.9328, 2.5733, 2.1493, 2.1295, 1.9906, 2.0082, 1.9875])

# 4000
times.append([26.5930, 15.5480, 13.2987, 9.6170, 9.0637, 8.4627, 8.3687, 7.9383])

# 8000
times.append([110.6343, 57.6750, 42.9657, 37.8757, 33.6563, 27.7037, 27.5733, 30.1933])

speedup      = []
efficiency   = []

def main():

	if not speedup_effiency_plot:
		df = pd.read_csv("../build/Diffuse.txt", sep=" ")
		df.sort_values(by=['t', 'x'], inplace=True)

		t = np.array(df.iloc[:, 0].tolist())
		x = np.array(df.iloc[:, 1].tolist())
		u = np.array(df.iloc[:, 2].tolist())

		fig = plt.figure(figsize=[12, 8]) 
		ax = fig.add_subplot(111, projection ='3d')

		im = ax.scatter(x, t, u, cmap='coolwarm', c=u)
		cbar = fig.colorbar(im, orientation='vertical')

		ax.set_xlabel(r"$x$", fontsize=18)
		ax.set_ylabel(r"$t$", fontsize=18)
		ax.set_zlabel(r"$u$", fontsize=18)
		ax.set_xlim([0, 1])
		ax.set_ylim([0, 1])
		ax.set_zlim([-1, 1.5])

		plt.show()

	else:
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

			ax[0].plot(processes, new_speedup,    label=f"n = {mesh_ns[i]}")
			ax[1].plot(processes, new_efficiency, label=f"n = {mesh_ns[i]}")



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