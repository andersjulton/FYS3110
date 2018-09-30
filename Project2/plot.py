import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('TkAgg')

with open("wave_21.txt", "r") as infile:
    n = int(infile.readline())
    u1 = np.zeros(n)
    for i in range(n):
        u1[i] = infile.readline()

with open("wave_22.txt", "r") as infile:
    n = int(infile.readline())
    u2 = np.zeros(n)
    for i in range(n):
        u2[i] = infile.readline()

with open("wave_23.txt", "r") as infile:
    n = int(infile.readline())
    u3 = np.zeros(n)
    for i in range(n):
        u3[i] = infile.readline()

with open("wave_24.txt", "r") as infile:
    n = int(infile.readline())
    u4 = np.zeros(n)
    for i in range(n):
        u4[i] = infile.readline()

x = np.linspace(0, 1, n)

plt.figure()
plt.plot(x, u1, color='red')
plt.plot(x, u2, color='green')
plt.plot(x, u3, color='blue')
plt.plot(x, u4, color='orange')
plt.show()
