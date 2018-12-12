import numpy as np

E = np.array([-8, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 0, 0, 0, 0, -8])
M = np.array([4, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, -2, -2, -2, -2 , -4])

def Z(T):
    return 2*np.exp(-8*T) + 2*np.exp(8*T) + 12

def mE(E, pow, T):
    sum = 0
    for i in range(len(E)):
        sum += E[i]**(pow)*np.exp(-E[i]*T)
    return sum/Z(T)

def CV(T):
    fac1 = mE(E, 2, T)
    fac2 = mE(E, 1, T)**2
    return (fac1 - fac2)/T**2


def mM(E, M, pow, T):
    sum = 0
    for i in range(len(E)):
        sum += np.abs(M[i])**(pow)*np.exp(-E[i]*T)
    return sum/Z(T)

def susc(T):
    fac1 = mM(E, M, 2, T)
    fac2 = mM(E, M, 1, T)**2
    return (fac1 - fac2)/T


E1 = -32*np.sinh(8.0)/(4*np.cosh(8.0) + 12)
E2 = 256*np.cosh(8.0)/(4*np.cosh(8.0) + 12)

print(E2 - E1**2)


print("Mean energy = %5.8f" %mE(E, 1, 2.0))
print("CV = %5.8f" %CV(1.0))
print("Mean magnetization = %5.8f" %mM(E, M, 1, 2.0))
print("Magnetic susceptibility = %5.8f" %susc(1.0))
