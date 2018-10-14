import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
	with open(filename, "r") as infile:
		n = int(infile.readline().strip())
		v = np.zeros(n)
		for i in range(n):
			v[i] = float(infile.readline().strip())
	return v

pos_Ex = read_file("pos_Ex.txt")
pos_Ey = read_file("pos_Ey.txt")
pos_Ez = read_file("pos_Ez.txt")
pos_Jx = read_file("pos_Jx.txt")
pos_Jy = read_file("pos_Jy.txt")
pos_Jz = read_file("pos_Jz.txt")

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_zlim(-1, 1)
ax.plot(pos_Ex, pos_Ey, pos_Ez, label='Earth')
ax.plot(pos_Jx, pos_Jy, pos_Jz, label='Jupiter')
ax.legend()
plt.show()
