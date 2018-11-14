import numpy as np
import matplotlib.pyplot as plt

def relError(a, b):
    return np.abs(np.abs(a - b)/a)

def Z(temp):
    Z = 2*np.exp(8/temp) + 2*np.exp(-8/temp) + 12
    return Z

def eMeanA(temp):
    E = -(1.0/Z(temp))*(16*np.exp(8/temp) - 16*np.exp(-8/temp))
    return E

def eSqMeanA(temp):
    ESq = (1.0/Z(temp))*(128*np.exp(8/temp) + 128*np.exp(-8/temp))
    return ESq

def mabsMeanA(temp):
    M = (1.0/Z(temp))*(8*np.exp(8/temp) + 16)
    return M

def mSqMeanA(temp):
    MSq = (1.0/Z(temp))*(32*np.exp(8/temp) + 32)
    return MSq

def CVA(temp):
    CV = (1.0/(temp*temp))*(eSqMeanA(temp) - eMeanA(1.0)*eMeanA(1.0))
    return CV

def suscA(temp):
    susc = (1.0/(temp))*(mSqMeanA(temp) - mabsMeanA(temp)*mabsMeanA(temp))
    return susc


def read_file(filename, n, m):
	a = np.fromfile(filename)
	return a.reshape(m, n)

def get_dim(filename):
    with open(filename+"_INFO.txt", "r") as infile:
        line1 = infile.readline()
        n = int(line1.split()[0])
        m = float(line1.split()[1])
        temp = float(line1.split()[2])
        m = int(m)
    return n, m, temp

def ex_b():
    n, m, temp = get_dim("ex_b0")
    a0 = read_file("ex_b0", m, n)
    a1 = read_file("ex_b1", m, n)
    a2 = read_file("ex_b2", m, n)
    a3 = read_file("ex_b3", m, n)
    a4 = read_file("ex_b4", m, n)

    a = (a0 + a1 + a2 + a3 + a4)/5

    x = np.linspace(10, 10*m, m)
    plotlim = 1000

    plt.figure(figsize=(12,7))
    plt.subplot(221)
    plt.plot(x[0:plotlim], a[0][0:plotlim])
    plt.axhline(color = 'orange', y=eMeanA(temp))
    plt.ylabel(r'$\langle E\rangle$')
    plt.xlabel("Monte Carlo cycles")

    plt.subplot(222)
    plt.plot(x[0:plotlim], a[1][0:plotlim])
    plt.axhline(color = 'orange', y=mabsMeanA(temp))
    plt.ylabel(r'$\langle |M|\rangle$')
    plt.xlabel("Monte Carlo cycles")

    plt.subplot(223)
    plt.plot(x[0:plotlim], a[2][0:plotlim])
    plt.axhline(color = 'orange', y=CVA(temp))
    plt.ylabel(r'$C_V$')
    plt.xlabel("Monte Carlo cycles")

    plt.subplot(224)
    plt.plot(x[0:plotlim], a[3][0:plotlim])
    plt.axhline(color = 'orange', y=suscA(temp))
    plt.ylabel(r'$\chi$')
    plt.xlabel("Monte Carlo cycles")

    plt.figure()
    plt.semilogy(x, relError(eMeanA(1.0), a[0]), label=r'$\langle E\rangle$')
    plt.semilogy(x, relError(mabsMeanA(1.0), a[1]), label=r'$\langle |M|\rangle$')
    plt.semilogy(x, relError(CVA(1.0), a[2]), label=r'$C_V$')
    plt.semilogy(x, relError(suscA(1.0), a[3]), label=r'$\chi$')
    plt.xlabel('Monte Carlo cycles', fontsize=15)
    plt.ylabel('Relative error', fontsize=15)
    plt.legend(fontsize=12)
    plt.savefig('b_error.pdf')

    plt.show()

def ex_c():
    filename1 = "ex_b_1_unord"
    filename2 = "ex_b_1_ord"
    filename3 = "ex_b_2.4_unord"
    filename4 = "ex_b_2.4_ord"
    n, m, temp = get_dim(filename1)
    values1 = read_file(filename1, m, n)
    values2 = read_file(filename2, m, n)
    values3 = read_file(filename3, m, n)
    values4 = read_file(filename4, m, n)
    x = np.linspace(10, 10*m, m)

    plt.figure(figsize=(12,10))
    plt.subplot(221)
    plt.plot(x, values1[0], label='Random start config')
    plt.plot(x, values2[0], label='Ordered start config')
    plt.xlabel("Monte Carlo cycles", fontsize=12)
    plt.ylabel(r'$\langle E\rangle$', fontsize=12)
    plt.title(r'T = 1.0')

    plt.subplot(222)
    plt.plot(x, values3[0], label='Random start config')
    plt.plot(x, values4[0], label='Ordered start config')
    plt.xlabel("Monte Carlo cycles", fontsize=12)
    plt.ylabel(r'$\langle E\rangle$', fontsize=12)
    plt.title(r'T = 2.4')

    plt.subplot(223)
    plt.plot(x, values1[3], label='Random start config')
    plt.plot(x, values2[3], label='Ordered start config')
    plt.xlabel("Monte Carlo cycles", fontsize=12)
    plt.ylabel(r'$\langle |M|\rangle$', fontsize=12)
    plt.title(r'T = 1.0')

    plt.subplot(224)
    plt.plot(x, values3[3])
    plt.plot(x, values4[3])
    plt.xlabel("Monte Carlo cycles", fontsize=12)
    plt.ylabel(r'$\langle |M|\rangle$', fontsize=12)
    plt.title(r'T = 2.4')


    plt.figure()
    plt.subplot(211)
    plt.plot(x, values1[4])
    plt.plot(x, values2[4])
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel("Accepted configurations")
    plt.title(r'T = 1.0')

    plt.subplot(212)
    plt.plot(x, values3[4])
    plt.plot(x, values4[4])
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel("Accepted configurations")
    plt.title(r'T = 2.4')

    plt.show()

def ex_d():
    filename1 = "ex_d_1"
    filename2 = "ex_d_2.4"

    n1, m1, temp1 = get_dim(filename1)
    values1 = read_file(filename1, m1, n1)
    n2, m2, temp2 = get_dim(filename2)
    values2 = read_file(filename2, m2, n2)

    sortindex1 = values1[0].argsort()
    sortindex2 = values2[0].argsort()

    print(values1[0][sortindex1])

    plt.figure()
    plt.subplot(121)
    plt.plot(values1[0][sortindex1], values1[1][sortindex1])
    plt.subplot(122)
    plt.plot(values2[0][sortindex2], values2[1][sortindex2])

    plt.show()

def ex_e():
    filename40 = "ex_e_40x40_0.005"
    filename60 = "ex_e_60x60_0.005"
    filename80 = "ex_e_80x80_0.005"
    filename100 = "ex_e_100x100_0.005"
    values40 = read_file(filename40, 21, 7)
    values60 = read_file(filename60, 21, 7)
    values80 = read_file(filename80, 21, 7)
    values100 = read_file(filename100, 21, 7)

    values = [values40, values60, values80, values100]

    mMax = np.zeros(4)
    suscMax = mMax.copy()
    CVMax = mMax.copy()
    x = np.linspace(2.2, 2.3, 21)

    for i in range(4):
        values[i][5] /= x**2
        values[i][6] /= x

    for i in range(4):
        mMax[i] = values[i][4][0]
        suscMax[i] = values[i][6].max()
        CVMax[i] = values[i][5].max()
        plt.plot(x, values[i][0])


    lsize = [40, 60, 80, 100]

    plt.show()

def plot_meanE():
    plt.plot(x, values40[0], label='40x40')
    plt.plot(x, values60[0], label='60x60')
    plt.plot(x, values80[0], label='80x80')
    plt.plot(x, values100[0], label='100x100')
    plt.legend()
    plt.show()

def plot_meanM():
    plt.plot(x, values40[4], label='40x40')
    plt.plot(x, values60[4], label='60x60')
    plt.plot(x, values80[4], label='80x80')
    plt.plot(x, values100[4], label='100x100')
    plt.legend()
    plt.show()

def plot_CV():
    plt.plot(x, values40[5]/x**2, label='40x40')
    plt.plot(x, values60[5]/x**2, label='60x60')
    plt.plot(x, values80[5]/x**2, label='80x80')
    plt.plot(x, values100[5]/x**2, label='100x100')
    plt.legend()
    plt.show()

def plot_susc():
    plt.plot(x, values40[6]/x, label='40x40')
    plt.plot(x, values60[6]/x, label='60x60')
    plt.plot(x, values80[6]/x, label='80x80')
    plt.plot(x, values100[6]/x, label='100x100')
    plt.legend()
    plt.show()

#plot_meanE()
#plot_meanM()
#plot_CV()
#plot_susc()
