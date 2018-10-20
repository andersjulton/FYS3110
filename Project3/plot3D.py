import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def readFileMatrix(filename):
	with open(filename, "r") as infile:
		line_1 = infile.readline()
		n = int(line_1.split()[0])
		m = int(line_1.split()[1])
		v = np.zeros((n, m))
		for i in range(n):
			line = infile.readline()
			for j in range(m):
				v[i][j] = float(line.split()[j])
	return v

def readFileArray(filename):
	with open(filename, "r") as infile:
		n = int(infile.readline().strip())
		v = np.zeros(n)
		for i in range(n):
			v[i] = float(infile.readline().strip())
	return v

def plotMany():
	names = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	for i, n in enumerate(names):
		ax.plot(pos_x[:,i], pos_y[:,i], pos_z[:,i], label = n, linewidth=0.5)
	ax.legend()
	plt.show()

def plotOne(planet):
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(pos_x, pos_y, pos_z, label = planet, linewidth=0.5)
	ax.legend()
	plt.show()

pos_x = readFileArray("ES_pos_x.txt")
pos_y = readFileArray("ES_pos_y.txt")


plt.plot(pos_x[0:100], pos_y[0:100])
plt.show()
