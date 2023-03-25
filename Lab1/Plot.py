#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


def analytic1(x, t):
	u = x*t - t *t / 2
	u += ((2*t - x)**2) / 8 + np.exp(-(t - x/2))

	return u

def analytic2(x, t):
	u = x*t - t *t / 2
	u += np.cos(np.pi * (2*t - x))
	
	return u


def main():
	df = pd.read_csv("../build/HeatEq.txt", sep=" ")
	df.sort_values(by=['t', 'x'], inplace=True)

	t = np.array(df.iloc[:, 0].tolist())
	x = np.array(df.iloc[:, 1].tolist())
	u = np.array(df.iloc[:, 2].tolist())

	
	df_par = pd.read_csv("../build/HeatEqTest.txt", sep=" ")
	df_par.sort_values(by=['t', 'x'], inplace=True)

	t_par = np.array(df_par.iloc[:, 0].tolist())
	x_par = np.array(df_par.iloc[:, 1].tolist())
	u_par = np.array(df_par.iloc[:, 2].tolist())


	fig = plt.figure(figsize=[12, 8]) 
	ax = fig.add_subplot(111, projection ='3d')

	t_an = np.linspace(0, 1, 100)
	x_an = np.linspace(0, 1, 100)

	#X, Y = np.meshgrid(x_an, t_an)
	#Z = analytic1(X, Y)
	#ax.contour3D(X, Y, Z, 50, cmap='binary')

	#Z = analytic2(X, Y)
	#ax.contour3D(X, Y, Z, 50, cmap='binary')


	im = ax.scatter(x_par, t_par, u_par, cmap='coolwarm', c=u_par)
	cbar = fig.colorbar(im, orientation='vertical')
	cbar.set_label('Temperature')
	#ax.plot(x_par, t_par, u - u_par)

	ax.set_xlabel(r"$x$", fontsize=18)
	ax.set_ylabel(r"$t$", fontsize=18)
	ax.set_zlabel(r"$u$", fontsize=18)
	ax.set_xlim([0, 1])
	ax.set_ylim([0, 1])
	ax.set_zlim([-1, 1.5])


	plt.show()



if __name__ == '__main__':
    main()