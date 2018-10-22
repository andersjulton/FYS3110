import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def get_dim(filename):
	with open(filename+"_DIM.txt", "r") as infile:
		line1 = infile.readline()
		n = int(line1.split()[0])
		m = float(line1.split()[1])
		m = int(m)
	return n, m

def read_file(filename, n, m):
	a = np.fromfile(filename)
	return a.reshape(m, n)

def plot_3D(names, filename):
	n, m = get_dim(filename)
	pos_x = read_file(filename + "_x", n, m)
	pos_y = read_file(filename + "_y", n, m)
	pos_z = read_file(filename + "_z", n, m)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	for i, n in enumerate(names):
		ax.plot(pos_x[i], pos_y[i], pos_z[i], label = n, linewidth=0.5)
	ax.legend()
	#plt.savefig(filename + ".png")
	#plt.close()
	plt.show()

def plot_2D_compare(filenames, figname, labels):
	for filename in filenames:
		n, m = get_dim(filename)
		pos_x = read_file(filename + "_x", n, m)
		pos_y = read_file(filename + "_y", n, m)
		plt.plot(pos_x, pos_y)
	#plt.legend(labels, fontsize=15)
	#plt.axis('equal')
	plt.xticks(fontsize=15)
	plt.yticks(fontsize=15)
	plt.grid()

def plotComp():
	plt.figure(figsize=(10,10))
	plt.subplot(221)
	plt.title(r'Comparison for $\Delta t = 1e-3$', fontsize=15)
	plt.ylabel(r'$y \ \ [AU]$', fontsize=15)
	plot_2D_compare(["earth_sun_euler_100", "earth_sun_verlet_100"], "earth_sun", ["Eulers", "Verlet"])
	plt.subplot(222)
	plt.title(r'Comparison for $\Delta t = 1e-4$', fontsize=15)
	plot_2D_compare(["earth_sun_euler_1000", "earth_sun_verlet_1000"], "earth_sun", ["Eulers", "Verlet"])
	plt.subplot(223)
	plt.title(r'Comparison for $\Delta t = 1e-5$', fontsize=15)
	plt.xlabel(r'$x \ \ [AU]$', fontsize=15)
	plt.ylabel(r'$y \ \ [AU]$', fontsize=15)
	plot_2D_compare(["earth_sun_euler_10000", "earth_sun_verlet_10000"], "earth_sun", ["Eulers", "Verlet"])
	plt.subplot(224)
	plt.title(r'Comparison for $\Delta t = 1e-6$', fontsize=15)
	plt.xlabel(r'$x  \ \ [AU]$', fontsize=15)
	plot_2D_compare(["earth_sun_euler_100000", "earth_sun_verlet_100000"], "earth_sun", ["Eulers", "Verlet"])
	plt.legend(["Eulers", "Verlet"], fontsize=15)
	plt.savefig("Comparison.pdf")
	plt.show()


names = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
plot_3D(names, "whole")
n, m = get_dim("whole")
