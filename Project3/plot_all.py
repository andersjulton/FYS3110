import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
	with open(filename, "r") as infile:
		line_1 = infile.readline()
		n = int(line_1.split()[0])
		m = int(line_1.split()[1])
		v = np.zeros((n, m))
		for i in range(n):
			line = infile.readline().split()
			for j in range(m):
				v[i][j] = float(line[j])
	return v

def plot_3D(names, filename):
	pos_x = read_file(filename + "_x.txt")
	pos_y = read_file(filename + "_y.txt")
	pos_z = read_file(filename + "_z.txt")
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	for i, n in enumerate(names):
		ax.plot(pos_x[:,i], pos_y[:,i], pos_z[:,i], label = names[i], linewidth=0.5)
	ax.legend()
	plt.savefig(filename + ".png")
	plt.close()

def plot_2D_compare(filenames, figname, labels, title):
	fig = plt.figure(1, figsize = (8,8))
	ax = fig.add_subplot(1, 1, 1)
	for filename in filenames:
		pos_x = read_file(filename + "_x.txt")
		pos_y = read_file(filename + "_y.txt")
		pos_z = read_file(filename + "_z.txt")
		ax.plot(pos_x[:,0], pos_y[:,0])
	plt.legend(labels, fontsize=15)
	plt.xlabel(r"$x [AU]$")
	plt.ylabel(r"$y [AU]$")
	plt.xticks(fontsize=15)
	plt.title(title)
	plt.yticks(fontsize=15)
	#plt.ylim(-50, 150)
	#plt.xlim(-1600, 50)
	plt.grid()
	plt.savefig(figname + ".png")
	plt.close()

plot_2D_compare(["Earth_Sun_min", "Earth_Sun_max"], "escapeVel", \
	[r"$v = 8.875 AU/yr$", r"$v = 8.89 AU/yr$"], "")