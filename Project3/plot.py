import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def read_file(filename):
	with open(filename, "r") as infile:
		n = int(infile.readline().strip())
		v = np.zeros(n)
		for i in range(n):
			v[i] = float(infile.readline().strip())
	return v

pos_x = read_file("pos_x.txt")
pos_y = read_file("pos_y.txt")

plt.plot(pos_x, pos_y, label='Earth')

pos_Ex = read_file("pos_Ex.txt")
pos_Ey = read_file("pos_Ey.txt")
pos_Jx = read_file("pos_Jx.txt")
pos_Jy = read_file("pos_Jy.txt")

plt.plot(pos_Ex, pos_Ey, label = 'Earth')
plt.plot(pos_Jx, pos_Jy, label = 'Jupiter')


plt.legend(loc = 'best')
plt.axis('equal')
plt.show()
