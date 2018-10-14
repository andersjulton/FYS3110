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

pos_Ex = read_file("Earth_x.txt")
pos_Ey = read_file("Earth_y.txt")
pos_Ez = read_file("Earth_z.txt")
pos_MAx = read_file("Mars_x.txt")
pos_MAy = read_file("Mars_y.txt")
pos_MAz = read_file("Mars_z.txt")
pos_Vx = read_file("Venus_x.txt")
pos_Vy = read_file("Venus_y.txt")
pos_Vz = read_file("Venus_z.txt")
pos_Px = read_file("Pluto_x.txt")
pos_Py = read_file("Pluto_y.txt")
pos_Pz = read_file("Pluto_z.txt")
pos_MEx = read_file("Mercury_x.txt")
pos_MEy = read_file("Mercury_y.txt")
pos_MEz = read_file("Mercury_z.txt")
pos_Jx = read_file("Jupiter_x.txt")
pos_Jy = read_file("Jupiter_y.txt")
pos_Jz = read_file("Jupiter_z.txt")

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_zlim(-0.5, 0.5)
ax.plot(pos_Ex, pos_Ey, pos_Ez, label='Earth')
ax.plot(pos_MAx, pos_MAy, pos_MAz, label='Mars')
ax.plot(pos_Vx, pos_Vy, pos_Vz, label='Venus')
ax.plot(pos_MEx, pos_MEy, pos_MEz, label='Mercury')
ax.plot(pos_Jx, pos_Jy, pos_Jz, label='Jupiter')
ax.legend()
plt.show()
