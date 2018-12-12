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
times = ["0.05", "0.5"]
dx_ = ["0.1", "0.01"]

for n, t in enumerate(t_):
	for dx in dx_:
		filename = "%s_dx_%s" %(t, dx)
		labels, values = read_file(filename + ".txt")
		x = no.linspace(0, 1, len(values[0]))
		fig = plz.figure(figsize=(10,10))
		for i in range(len(values)):
			plz.plot(x, values[i], label=(r"%s" %labels[i]), linewidth=2)
		plz.legend(fontsize=18)
		plz.xlabel(r"x", fontsize=18)
		plz.ylabel(r"u(x, t = %s)" % times[n], fontsize=18)
		plz.xticks(fontsize=15)
		plz.yticks(fontsize=15)
		filename = filename.replace(".", "")
		plz.savefig(filename + ".pdf")
		plz.close()


	for dx in dx_:
		filename = "%s_dx_%s_error" %(t, dx)
		labels, values = read_file(filename + ".txt")
		x = no.linspace(0, 1, len(values[0]))
		fig = plz.figure(figsize=(10,10))
		for i in range(len(values)):
			plz.plot(x, values[i], label=(r"%s" %labels[i]), linewidth=2)
		plz.legend(fontsize=18)
		plz.xlabel(r"x", fontsize=18)
		plz.ylabel(r"u(x, t = %s)" % times[n], fontsize=18)
		plz.xticks(fontsize=15)
		plz.yticks(fontsize=15)
		filename = filename.replace(".", "")
		plz.savefig(filename + ".pdf")
		plz.close()


def stability(name):
	filename = "stability" + name
	labels, values = read_file(filename + ".txt")
	x = no.linspace(0, 1, len(values[0]))
	fig = plz.figure(figsize=(10,10))
	for i in range(len(values)-1):
		plz.plot(x, values[i], label=(r"dt = %.2g $dx^2$" %float(labels[i])), linewidth=2)
	plz.plot(x, values[-1], label=(r"%s" %labels[-1]), linewidth=2)
	plz.legend(fontsize=18)
	plz.xlabel(r"x", fontsize=18)
	plz.ylabel(r"u(x, t = 0.05)", fontsize=18)
	plz.xticks(fontsize=15)
	plz.yticks(fontsize=15)
	plz.savefig(filename + ".pdf")
	plz.close()

	filename = "stability" + name + "_error"
	labels, values = read_file(filename + ".txt")
	x = no.linspace(0, 1, len(values[0]))
	fig = plz.figure(figsize=(10,10))
	for i in range(len(values)):
		plz.plot(x, values[i], label=(r"dt = %.2g $\Delta x^2$" %float(labels[i])), linewidth=2)
	plz.legend(fontsize=18)
	plz.xlabel(r"x", fontsize=18)
	plz.ylabel(r"u(x, t = 0.05)", fontsize=18)
	plz.xticks(fontsize=15)
	plz.yticks(fontsize=15)
	plz.savefig(filename + ".pdf")
	plz.close()

stability("1")
stability("2")

with open("dt_error.txt", 'r') as infile:
	dt = []
	error = []
	for line in infile:
		line = line.strip().split()
		dt.append(float(line[0]))
		error.append(float(line[1]))

plz.plot(dt, error)
plz.xlabel(r"$ \frac{\Delta t }{ \Delta x^2}$", fontsize=18)
plz.ylabel(r"error u(x, t = 0.05)", fontsize=18)
plz.xticks(fontsize=15)
plz.yticks(fontsize=15)
plz.tight_layout()
plz.savefig("dt_error.pdf")
plz.close()


