import matplotlib.pyplot as plt
import numpy as np

n = 100
A = [[] for i in range(n)]
for i in range(n):
    for j in range(n):
        A[i].append((i*n + j))
        
x = np.linspace(0, 1, n)


#cmap = plt.get_cmap('magma')

fig = plt.figure(1, figsize = (8,12))

ax1 = fig.add_subplot(2, 1, 1)
cmap = plt.get_cmap('inferno')
im = ax1.pcolor(A, cmap=cmap)
fig.colorbar(im, ax=ax1)
plt.xlabel(r"x")
plt.ylabel(r"y")

ax2 = fig.add_subplot(2, 1, 2)
cmap = plt.get_cmap('Greys')
im = ax2.pcolor(A, cmap=cmap)
fig.colorbar(im, ax=ax2)
plt.xlabel(r"x")
plt.ylabel(r"y")


fig.tight_layout()

plt.show()