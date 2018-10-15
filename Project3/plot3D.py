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
			line = infile.readline()
			for j in range(m):
				v[i][j] = float(line.split()[j])
	return v

pos_x = read_file("pos_x.txt")
pos_y = read_file("pos_y.txt")
pos_z = read_file("pos_z.txt")


names = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

fig = plt.figure()
ax = fig.gca(projection='3d')
for i, n in enumerate(names):
	ax.plot(pos_x[:,i], pos_y[:,i], pos_z[:,i], label = n, linewidth=0.5)
ax.legend()
plt.show()
