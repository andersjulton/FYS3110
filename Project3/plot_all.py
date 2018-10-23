import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def get_dim(filename):
	with open(filename+"_DIM.txt", "r") as infile:
		line1 = infile.readline()
		n = int(line1.split()[0])
		m = float(line1.split()[1])
		m = int(m)
	return n, m

def read_file(filename, n, m):
	a = np.fromfile(filename)
	return a.reshape(m, n)

def plot_3D(names, filename):
	n, m = get_dim(filename)
	pos_x = read_file(filename + "_x", n, m)
	pos_y = read_file(filename + "_y", n, m)
	pos_z = read_file(filename + "_z", n, m)
	maxX = np.max(pos_x)
	maxY = np.max(pos_y)
	max = np.max([maxX, maxY])
	max = 1.2
	fig = plt.figure(figsize=(10,10))
	ax = fig.gca(projection='3d')
	ax.set_xlim(-max, max)
	ax.set_ylim(-max, max)
	ax.set_zlim(-max, max)
	ax.set_xlabel(r'$x \ \ [AU]$', fontsize=15)
	ax.set_ylabel(r'$y \ \ [AU]$', fontsize=15)
	ax.set_zlabel(r'$z \ \ [AU]$', fontsize=15)
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False

	# Now set color to white (or whatever is "invisible")
	ax.xaxis.pane.set_edgecolor('w')
	ax.yaxis.pane.set_edgecolor('w')
	ax.zaxis.pane.set_edgecolor('w')

	for i, n in enumerate(names):
		ax.plot(pos_x[i], pos_y[i], pos_z[i], label = n, linewidth=0.5)
	ax.legend()
	plt.savefig(filename + "inner.pdf")
	#plt.close()
	plt.show()

def plot_SEJ(names, filename):
	n, m = get_dim(filename)
	r = 0.00464913034*20
	pos_x = read_file(filename + "_x", n, m)
	pos_y = read_file(filename + "_y", n, m)
	pos_z = read_file(filename + "_z", n, m)
	maxX = np.max(pos_x)
	maxY = np.max(pos_y)
	max = np.max([maxX, maxY])
	fig = plt.figure(figsize=(10,10))
	ax = fig.gca(projection='3d')
	ax.set_xlim(-max, max)
	ax.set_ylim(-max, max)
	ax.set_zlim(-max, max)
	ax.set_xlabel(r'$x \ \ [AU]$', fontsize=15)
	ax.set_ylabel(r'$y \ \ [AU]$', fontsize=15)
	ax.set_zlabel(r'$z \ \ [AU]$', fontsize=15)
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False

	for i, n in enumerate(names):
		ax.plot(pos_x[i], pos_y[i], pos_z[i], label = n, linewidth=0.5)
	ax.legend()
	#plt.savefig(filename + ".pdf")
	#plt.close()
	plt.show()

def plot_EJ_mass_compare(names, filenames, filename):
	r = 0.00464913034*20
	fig = plt.figure(figsize=(10,10))
	ax = fig.gca(projection='3d')
	ax.set_xlabel(r'$x \ \ [AU]$', fontsize=15)
	ax.set_ylabel(r'$y \ \ [AU]$', fontsize=15)
	ax.set_zlabel(r'$z \ \ [AU]$', fontsize=15)
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	u = np.linspace(0, 2 * np.pi, 50)
	v = np.linspace(0, np.pi, 50)
	x = r * np.outer(np.cos(u), np.sin(v))
	y = r * np.outer(np.sin(u), np.sin(v))
	z = r * np.outer(np.ones(np.size(u)), np.cos(v))
	for k in filenames:
		n, m = get_dim(k)
		pos_x = read_file(k + "_x", n, m)
		pos_y = read_file(k + "_y", n, m)
		pos_z = read_file(k + "_z", n, m)
		maxX = np.max(pos_x)
		maxY = np.max(pos_y)
		max = 5
		for i, n in enumerate(names):
			ax.plot(pos_x[i], pos_y[i], pos_z[i], label = n, linewidth=0.5)
	ax.legend()
	ax.set_xlim(-max, max)
	ax.set_ylim(-max, max)
	ax.set_zlim(-max, max)
	# Plot the surface
	ax.plot_surface(x, y, z, color='yellow')
	ax.scatter(0.0,0.0,0.0, color='yellow', label='Sun')
	# Now set color to white (or whatever is "invisible")
	ax.xaxis.pane.set_edgecolor('w')
	ax.yaxis.pane.set_edgecolor('w')
	ax.zaxis.pane.set_edgecolor('w')
	#plt.savefig(filename + ".pdf")
	#plt.close()
	plt.show()

def plot_EJ_compare(names, filenames, filename):
	r = 0.00464913034*20
	fig = plt.figure(figsize=(10,10))
	ax = fig.gca(projection='3d')
	ax.set_xlabel(r'$x \ \ [AU]$', fontsize=15)
	ax.set_ylabel(r'$y \ \ [AU]$', fontsize=15)
	ax.set_zlabel(r'$z \ \ [AU]$', fontsize=15)
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	u = np.linspace(0, 2 * np.pi, 50)
	v = np.linspace(0, np.pi, 50)
	x = r * np.outer(np.cos(u), np.sin(v))
	y = r * np.outer(np.sin(u), np.sin(v))
	z = r * np.outer(np.ones(np.size(u)), np.cos(v))
	for k in filenames:
		n, m = get_dim(k)
		pos_x = read_file(k + "_x", n, m)
		pos_y = read_file(k + "_y", n, m)
		pos_z = read_file(k + "_z", n, m)
		maxX = np.max(pos_x)
		maxY = np.max(pos_y)
		max = 5
		for i, n in enumerate(names):
			ax.plot(pos_x[i], pos_y[i], pos_z[i], label = n, linewidth=0.5)
	ax.legend()
	ax.set_xlim(-max, max)
	ax.set_ylim(-max, max)
	ax.set_zlim(-max, max)
	# Plot the surface
	ax.plot_surface(x, y, z, color='yellow')
	ax.scatter(0.0,0.0,0.0, color='yellow', label='Sun')
	# Now set color to white (or whatever is "invisible")
	ax.xaxis.pane.set_edgecolor('w')
	ax.yaxis.pane.set_edgecolor('w')
	ax.zaxis.pane.set_edgecolor('w')
	#plt.savefig(filename + ".pdf")
	#plt.close()
	plt.show()

def plot_2D_compare(filenames):
	for filename in filenames:
		n, m = get_dim(filename)
		pos_x = read_file(filename + "_x", n, m)
		pos_y = read_file(filename + "_y", n, m)
		plt.plot(pos_x, pos_y)
	#plt.legend(labels, fontsize=15)
	#plt.axis('equal')
	plt.xticks(fontsize=15)
	plt.yticks(fontsize=15)
	plt.grid()
	plt.show()

def plotComp():
	plt.figure(figsize=(10,10))
	plt.subplot(221)
	plt.title(r'Comparison for $\Delta t = 1e-3$', fontsize=15)
	plt.ylabel(r'$y \ \ [AU]$', fontsize=15)
	plot_2D_compare(["earth_sun_euler_100", "earth_sun_verlet_100"], "earth_sun", ["Eulers", "Verlet"])
	plt.subplot(222)
	plt.title(r'Comparison for $\Delta t = 1e-4$', fontsize=15)
	plot_2D_compare(["earth_sun_euler_1000", "earth_sun_verlet_1000"], "earth_sun", ["Eulers", "Verlet"])
	plt.subplot(223)
	plt.title(r'Comparison for $\Delta t = 1e-5$', fontsize=15)
	plt.xlabel(r'$x \ \ [AU]$', fontsize=15)
	plt.ylabel(r'$y \ \ [AU]$', fontsize=15)
	plot_2D_compare(["earth_sun_euler_10000", "earth_sun_verlet_10000"], "earth_sun", ["Eulers", "Verlet"])
	plt.subplot(224)
	plt.title(r'Comparison for $\Delta t = 1e-6$', fontsize=15)
	plt.xlabel(r'$x  \ \ [AU]$', fontsize=15)
	plot_2D_compare(["earth_sun_euler_100000", "earth_sun_verlet_100000"], "earth_sun", ["Eulers", "Verlet"])
	plt.legend(["Eulers", "Verlet"], fontsize=15)
	plt.savefig("Comparison.pdf")
	plt.show()

def plotCompBeta():
	labels = [r"v = 6.8 AU/yr", r"v = 7.4 AY/yr", r"v = 8.0 AU/yr", r"v = 8.8 AU/yr"]
	plt.figure(figsize=(12,12))
	plt.subplot(221)
	plt.title(r' $\beta = 2$', fontsize=15)
	plt.ylabel(r'$y \ \ [AU]$', fontsize=15)
	plot_2D_compare(["Beta_0_1", "Beta_0_2", "Beta_0_3", "Beta_0_4"], "", labels)
	plt.xlim(-20, 5)
	plt.ylim(-4, 8)
	plt.subplot(222)
	plt.title(r'$\beta = 2.33$', fontsize=15)
	plot_2D_compare(["Beta_1_1", "Beta_1_2", "Beta_1_3", "Beta_1_4"], "", labels)
	plt.xlim(-40, 10)
	plt.ylim(-10, 40)
	plt.subplot(223)
	plt.title(r'$\beta = 2.67$', fontsize=15)
	plt.xlabel(r'$x \ \ [AU]$', fontsize=15)
	plt.ylabel(r'$y \ \ [AU]$', fontsize=15)
	plot_2D_compare(["Beta_2_1", "Beta_2_2", "Beta_2_3", "Beta_2_4"], "", labels)
	plt.xlim(-50, 10)
	plt.ylim(-10, 60)
	plt.subplot(224)
	plt.title(r'$\beta = 3$', fontsize=15)
	plt.xlabel(r'$x  \ \ [AU]$', fontsize=15)
	plot_2D_compare(["Beta_3_1", "Beta_3_2", "Beta_3_3", "Beta_3_4"], "", labels)
	plt.xlim(-50, 5)
	plt.ylim(-40, 60)
	plt.legend(labels, fontsize=15)
	plt.savefig("ComparisonBeta.pdf")
	plt.close()

def plot_E_compare():
	names = ["Earth", "Jupiter"]
	filename = "earth_jupiter1x"
	n, m = get_dim(filename)
	nE, mE = get_dim("earth_sun")
	r = 0.00464913034*5
	pos_x = read_file(filename + "_x", n, m)
	pos_y = read_file(filename + "_y", n, m)
	pos_z = read_file(filename + "_z", n, m)
	Ex = read_file("earth_sun_x", nE, mE)
	Ey = read_file("earth_sun_y", nE, mE)
	Ez = read_file("earth_sun_z", nE, mE)
	maxX = np.max(Ex)
	maxY = np.max(Ey)
	max = np.max([maxX, maxY])*0.7
	fig = plt.figure(figsize=(10,10))
	ax = fig.gca(projection='3d')
	ax.set_xlim(-max, max)
	ax.set_ylim(-max, max)
	ax.set_zlim(-max, max)
	ax.set_xlabel(r'$x \ \ [AU]$', fontsize=15)
	ax.set_ylabel(r'$y \ \ [AU]$', fontsize=15)
	ax.set_zlabel(r'$z \ \ [AU]$', fontsize=15)
	ax.xaxis.pane.fill = False
	ax.yaxis.pane.fill = False
	ax.zaxis.pane.fill = False
	u = np.linspace(0, 2 * np.pi, 50)
	v = np.linspace(0, np.pi, 50)
	x = r * np.outer(np.cos(u), np.sin(v))
	y = r * np.outer(np.sin(u), np.sin(v))
	z = r * np.outer(np.ones(np.size(u)), np.cos(v))

	ax.plot_surface(x, y, z, color='yellow')
	ax.scatter(0.0,0.0,0.0, color='yellow', label='Sun')

	ax.xaxis.pane.set_edgecolor('w')
	ax.yaxis.pane.set_edgecolor('w')
	ax.zaxis.pane.set_edgecolor('w')

	ax.plot(pos_x[0], pos_y[0], pos_z[0], label = "Earth w/Jupiter", linewidth=0.5)
	ax.plot(Ex[0], Ey[0], Ez[0], label = "Earth wo/Jupiter", linewidth=0.5)
	ax.legend()
	#plt.savefig('earthJupiterCompare.pdf')
	plt.show()

def plot_energyAngular_whole():
	E = np.fromfile("whole_energy")
	L = np.fromfile("whole_angular")
	x = np.linspace(0, 1000, 1000)
	plt.figure(figsize=(15,5))
	plt.subplot(121)
	plt.xticks(fontsize=15)
	plt.yticks(fontsize=15)
	plt.xlabel('Years', fontsize=15)
	plt.ylabel(r'$E/E_0$', fontsize=15)
	plt.plot(x, E/E[0])
	plt.title('Change in total energy over time', fontsize=13)
	plt.subplot(122)
	plt.plot(x, L/L[0])
	plt.xticks(fontsize=15)
	plt.yticks(fontsize=15)
	plt.xlabel('Years', fontsize=15)
	plt.ylabel(r'$L/L_0$', fontsize=15)
	plt.title('Change in total angular momentum over time', fontsize=13, loc='right')
	plt.savefig('whole_E_L.pdf')
	plt.show()
names = ["Sun", "Mercury", "Venus", "Earth", "Mars"]#, "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
#plot_3D(names, "whole")
#plot_E_compare()
#plot_EJ_mass_compare(["Earth", "Jupiter"], ['earth_jupiter1x','earth_jupiter10x', 'earth_jupiter1000x'], "earth_jupiter_compare.pdf")
#plot_SEJ(["Sun", "Earth", "Jupiter"], "sunEarthJupiter")
#plot_2D_compare(["earthSun", "earth_jupiter1x"])


'''
filenameE = "earthSun"
filenameEJ = "earthJupiter"

nE, mE = get_dim(filenameE)
nEJ, mEJ = get_dim(filenameEJ)
Epos_x = read_file(filenameE + "_x", nE, mE)
Epos_y = read_file(filenameE + "_y", nE, mE)
Epos_z = read_file(filenameE + "_z", nE, mE)
EJpos_x = read_file(filenameEJ + "_x", nEJ, mEJ)
EJpos_y = read_file(filenameEJ + "_y", nEJ, mEJ)
EJpos_z = read_file(filenameEJ + "_z", nEJ, mEJ)

xdiff = np.abs(Epos_x - EJpos_x[0])
ydiff = np.abs(Epos_y - EJpos_y[0])
zdiff = np.abs(Epos_z - EJpos_z[0])

x = np.linspace(0, 100, nE)
plt.plot(x[0:50*nE/100], xdiff[0][0:50*nE/100], label= 'xdiff')
#plt.plot(ydiff[0], label= 'ydiff')
#plt.plot(zdiff[0], label= 'zdiff')
plt.xlabel(r'Time [Years]')
plt.show()'''

E = np.fromfile("earthJupiter_energy")
L = np.fromfile("earthJupiter_angular")
x = np.linspace(0, 100, 1000)
plt.figure(figsize=(15,5))
plt.subplot(121)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Years', fontsize=15)
plt.ylabel(r'$E/E_0$', fontsize=15)
plt.plot(x, E/E[0])
plt.title('Change in total energy over time', fontsize=13)
plt.subplot(122)
plt.plot(x, L/L[0])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Years', fontsize=15)
plt.ylabel(r'$L/L_0$', fontsize=15)
plt.title('Change in total angular momentum over time', fontsize=13, loc='right')
#plt.savefig('EJ_E_L.pdf')
plt.show()
