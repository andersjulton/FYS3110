import numpy as np
import matplotlib.pyplot as plt

with open("n_10.txt", "r") as infile:
    n_10 = int(infile.readline())
    v_10 = np.zeros(n_10)
    for i in range(n_10):
        v_10[i] = infile.readline()

with open("n_100.txt", "r") as infile:
    n_100 = int(infile.readline())
    v_100 = np.zeros(n_100)
    for i in range(n_100):
        v_100[i] = infile.readline()

with open("n_1000.txt", "r") as infile:
    n_1000 = int(infile.readline())
    v_1000 = np.zeros(n_1000)
    for i in range(n_1000):
        v_1000[i] = infile.readline()


x_10 = np.linspace(0, 1, n_10)
x_100 = np.linspace(0, 1, n_100)
x_1000 = np.linspace(0, 1, n_1000)
u = 1 - (1 - np.exp(-10))*x_1000 - np.exp(-10*x_1000)

plt.plot(x_1000, u, label="Exact solution")
plt.plot(x_10, v_10, label="n = 10")
plt.plot(x_100, v_100, label="n = 100")
plt.plot(x_1000, v_1000, label="n = 1000")

plt.legend(loc="best", fontsize=12)
plt.xlabel("x", fontsize=12)
plt.ylabel("u(x)", fontsize=12)
plt.savefig("compare.png")
plt.show()
