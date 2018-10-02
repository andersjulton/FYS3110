import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('TkAgg')

#Python script for plotting wave functions

with open("wave_1.txt", "r") as infile:
    n = int(infile.readline())
    u1 = np.zeros(n)
    for i in range(n):
        u1[i] = infile.readline()

with open("wave_2.txt", "r") as infile:
    n = int(infile.readline())
    u2 = np.zeros(n)
    for i in range(n):
        u2[i] = infile.readline()

with open("wave_3.txt", "r") as infile:
    n = int(infile.readline())
    u3 = np.zeros(n)
    for i in range(n):
        u3[i] = infile.readline()

with open("wave_4.txt", "r") as infile:
    n = int(infile.readline())
    u4 = np.zeros(n)
    for i in range(n):
        u4[i] = infile.readline()
with open("wave_1n.txt", "r") as infile:
    n = int(infile.readline())
    u1n = np.zeros(n)
    for i in range(n):
        u1n[i] = infile.readline()

with open("wave_2n.txt", "r") as infile:
    n = int(infile.readline())
    u2n = np.zeros(n)
    for i in range(n):
        u2n[i] = infile.readline()

with open("wave_3n.txt", "r") as infile:
    n = int(infile.readline())
    u3n = np.zeros(n)
    for i in range(n):
        u3n[i] = infile.readline()

with open("wave_4n.txt", "r") as infile:
    n = int(infile.readline())
    u4n = np.zeros(n)
    for i in range(n):
        u4n[i] = infile.readline()

x = np.linspace(0, 5, n)

plt.figure(figsize=(20,15))
plt.subplot(121)
plt.plot(x, u1, color='red', label=r'$\omega_r = 0.01$')
plt.plot(x, u2, color='green', label=r'$\omega_r = 0.5$')
plt.plot(x, u3, color='blue', label=r'$\omega_r = 1$')
plt.plot(x, u4, color='orange', label=r'$\omega_r = 5$')
plt.ylabel(r'$\psi_0(\rho)$')
plt.xlabel(r'$\rho$')
plt.legend(loc='best')
plt.subplot(122)
plt.plot(x, u1n, color='red', label=r'$\omega_r = 0.01$')
plt.plot(x, u2n, color='green', label=r'$\omega_r = 0.5$')
plt.plot(x, u3n, color='blue', label=r'$\omega_r = 1$')
plt.plot(x, u4n, color='orange', label=r'$\omega_r = 5$')
plt.ylabel(r'$\psi_0(\rho)$')
plt.xlabel(r'$\rho$')
plt.legend(loc='best')
#plt.savefig('two_electrons.png')
plt.show()
