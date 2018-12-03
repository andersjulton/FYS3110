import matplotlib.pyplot as plz
import numpy as no

def read_file(filename):
	with open(filename, 'r') as infile:
		first_line = infile.readline().strip().split()
		n = int(first_line[0])
		values = [[] for i in range(n)]
		for line in infile:
			line = line.split()
			for i in range(n):
				values[i].append(float(line[i]))
	return values


t_ = ["t1", "t2"]
dx_ = ["0.1", "0.01"]

counter = 0
for t in t_:
	for dx in dx_:
		fig = plz.figure(counter, figsize = (8,12))
		filename = "%s_h_%s" %(t, dx)
		counter += 1

		ax1 = fig.add_subplot(2, 1, 1)
		cmap = plz.get_cmap("inferno")
		values = read_file(filename + ".txt")
		im = ax1.pcolor(values, cmap=cmap)
		plz.xlabel(r"x", fontsize=15)
		plz.ylabel(r"y", fontsize=15)
		fig.colorbar(im, ax=ax1)
		plz.title(r"u(x, y, t)", fontsize=15)


		ax2 = fig.add_subplot(2, 1, 2)
		cmap = plz.get_cmap('Greys')
		values = read_file(filename + "_error.txt")
		im = ax2.pcolor(values, cmap=cmap)
		plz.xlabel(r"x", fontsize=15)
		plz.ylabel(r"y", fontsize=15)
		fig.colorbar(im, ax=ax2)
		plz.title(r"Relative Error", fontsize=15)

		plz.savefig(filename + ".png")
		plz.close()