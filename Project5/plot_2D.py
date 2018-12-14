import matplotlib.pyplot as plz
import numpy as no

def read_file(filename):
	with open(filename, 'r') as infile:
		first_line = infile.readline().strip().split()
		values = []
		for line in infile:
			line = [float(value) for value in line.split()]
			values.append(line)
	return values, int(first_line[0])


t_ = ["t1", "t2", "t3"]
times = ["0.005", "0.05", "0.5"]
dx_ = ["0.1", "0.01"]
counter = 0
for n, t in enumerate(t_):
	for dx in dx_:
		counter += 1
		fig = plz.figure(counter, figsize = (12,10))
		filename = "%s_h_%s" %(t, dx)
		ax1 = fig.add_subplot(1, 1, 1)
		cmap = plz.get_cmap("inferno")
		values, dim = read_file(filename + ".txt")
		X, Y = no.meshgrid(no.linspace(0, 1, dim+1), no.linspace(0, 1, dim+1))
		im = ax1.pcolormesh(X, Y, values, cmap=cmap, vmin=0, vmax=1)
		plz.xlabel(r"x", fontsize=20)
		plz.ylabel(r"y", fontsize=20)
		plz.xticks(fontsize=18)
		plz.yticks(fontsize=18)
		fig.colorbar(im, ax=ax1)
		plz.title(r"u(x, y, t = %s)" % times[n], fontsize=18)
		filename = filename.replace(".", "")
		plz.savefig(filename + ".pdf")
		plz.close()

		counter += 1
		filename = "%s_h_%s_error" %(t, dx)
		fig = plz.figure(counter, figsize = (12,10))
		ax2 = fig.add_subplot(1, 1, 1)
		cmap = plz.get_cmap('Greys')
		values, dim = read_file(filename + ".txt")
		im = ax2.pcolormesh(X, Y, values, cmap=cmap)
		plz.xlabel(r"x", fontsize=18)
		plz.ylabel(r"y", fontsize=18)
		plz.xticks(fontsize=15)
		plz.yticks(fontsize=15)
		fig.colorbar(im, ax=ax2)

		filename = filename.replace(".", "")
		plz.savefig(filename + ".pdf")
		plz.close()
