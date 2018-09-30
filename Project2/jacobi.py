import matplotlib.pyplot as plt 
import numpy as np

def read_file(filename):
	with open(filename, "r") as infile:
		n = int(infile.readline().strip())
		v = np.zeros(n)
		for i in range(n):
			v[i] = int(infile.readline().strip())
	return v

n = read_file("n.txt")
plt.figure(figsize=(12,8))
plt.plot(n, read_file("iter.txt"), ".", linewidth=3)
plt.xlabel("n", fontsize=20)
plt.ylabel("iterations", fontsize=20)
plt.ticklabel_format(style='sci', axis='y',  scilimits=(0,0))
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.savefig("iter.pdf")
plt.close()