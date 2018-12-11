import matplotlib.pyplot as plz
import numpy as no

def read_file(filename):
	with open(filename, 'r') as infile:
		first_line = infile.readline().strip().split()
		values = []
		for line in infile:
			line = [float(value) for value in line.split()]
			values.append(line)
	return values


t_ = ["t1", "t2", "t3"]
times = ["0.05", "0.05", "0.5"]
dx_ = ["0.1", "0.01"]
counter = 0
for n, t in enumerate(t_):
	for dx in dx_:
		counter += 1
		fig = plz.figure(counter, figsize = (12,10))
		filename = "%s_h_%s" %(t, dx)
		ax1 = fig.add_subplot(1, 1, 1)
		cmap = plz.get_cmap("inferno")
		values = read_file(filename + ".txt")
		im = ax1.pcolor(values, cmap=cmap, vmin=0, vmax=1)
		plz.xlabel(r"x", fontsize=15)
		plz.ylabel(r"y", fontsize=15)
		fig.colorbar(im, ax=ax1)
		plz.title(r"u(x, y, t = %s)" % times[n], fontsize=15)

		filename = filename.replace(".", "")
		plz.savefig(filename + ".png")
		plz.close()

		counter += 1
		filename = "%s_h_%s_error" %(t, dx)
		fig = plz.figure(counter, figsize = (12,10))
		ax2 = fig.add_subplot(1, 1, 1)
		cmap = plz.get_cmap('Greys')
		values = read_file(filename + ".txt")
		im = ax2.pcolor(values, cmap=cmap)
		plz.xlabel(r"x", fontsize=15)
		plz.ylabel(r"y", fontsize=15)
		fig.colorbar(im, ax=ax2)

		filename = filename.replace(".", "")
		plz.savefig(filename + ".png")
		plz.close()
