import matplotlib.pyplot as plz
import numpy as no

def read_file(filename):
	with open(filename, 'r') as infile:
		labels = infile.readline().strip().split("\t")
		n = len(labels)
		values = [[] for i in range(n)]
		for line in infile:
			line = line.split()
			for i in range(n):
				values[i].append(float(line[i]))
	return labels, values


t_ = ["t1", "t2"]
times = [0.05, 0.5]
dx_ = ["0.1", "0.01"]

for n, t in enumerate(t_):
	for dx in dx_:
		filename = "%s_dx_%s" %(t, dx)
		labels, values = read_file(filename + ".txt")
		x = no.linspace(0, 1, len(values[0]))
		fig = plz.figure(figsize=(10,10))
		for i in range(len(values)):
			plz.plot(x, values[i], label=(r"%s" %labels[i]), linewidth=2)
		plz.legend(fontsize=15)
		plz.xlabel(r"x", fontsize=15)
		plz.ylabel(r"u(x, t = %1.2f)" % times[n], fontsize=15)
		filename = filename.replace(".", "")
		plz.savefig(filename + ".png")
		plz.close()


	for dx in dx_:
		filename = "%s_dx_%s_error" %(t, dx)
		labels, values = read_file(filename + ".txt")
		x = no.linspace(0, 1, len(values[0]))
		fig = plz.figure(figsize=(10,10))
		for i in range(len(values)):
			plz.plot(x, values[i], label=(r"%s" %labels[i]), linewidth=2)
		plz.legend(fontsize=15)
		plz.xlabel(r"x", fontsize=15)
		plz.ylabel(r"u(x, t = %1.2f)" % times[n], fontsize=15)
		filename = filename.replace(".", "")
		plz.savefig(filename + ".png")
		plz.close()
