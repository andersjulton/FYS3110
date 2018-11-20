import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def relError(a, b):
    return np.abs(np.abs(a - b)/a)

def Z(temp):
    Z = 2*np.exp(8/temp) + 2*np.exp(-8/temp) + 12
    return Z

def Emean(temp):
    E = -(1.0/Z(temp))*(16*np.exp(8/temp) - 16*np.exp(-8/temp))
    return E

def E2Mean(temp):
    ESq = (1.0/Z(temp))*(128*np.exp(8/temp) + 128*np.exp(-8/temp))
    return ESq

def MAmean(temp):
    M = (1.0/Z(temp))*(8*np.exp(8/temp) + 16)
    return M

def M2mean(temp):
    MSq = (1.0/Z(temp))*(32*np.exp(8/temp) + 32)
    return MSq

def CVA(temp):
    CV = (1.0/(temp*temp))*(E2Mean(temp) - Emean(temp)*Emean(temp))
    return CV

def suscA(temp):
    susc = (1.0/(temp))*(M2mean(temp) - MAmean(temp)*MAmean(temp))
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

    plt.figure(figsize=(8,6))
    plt.loglog(x, relError(Emean(1.0), a[0]), label=r'$\langle E\rangle$')
    plt.loglog(x, relError(MAmean(1.0), a[1]), label=r'$\langle |M|\rangle$')
    plt.loglog(x, relError(CVA(1.0), a[2]), label=r'$C_V$')
    plt.loglog(x, relError(suscA(1.0), a[3]), label=r'$\chi$')
    plt.xlabel('Monte Carlo cycles', fontsize=15)
    plt.ylabel('Relative error', fontsize=15)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=14)
    plt.grid(True)
    plt.savefig('b_error.pdf')
    plt.show()

def ex_c():
    filename1 = "ex_c_1_unord"
    filename2 = "ex_c_1_ord"
    filename3 = "ex_c_2.4_unord"
    filename4 = "ex_c_2.4_ord"
    n, m, temp = get_dim(filename1)
    values1 = read_file(filename1, m, n)
    values2 = read_file(filename2, m, n)
    values3 = read_file(filename3, m, n)
    values4 = read_file(filename4, m, n)
    x = np.linspace(10, 10*m, m)
    plotlim = 500

    plt.figure(figsize=(12,11))
    plt.subplot(221)
    plt.title(r'T = 1.0')
    plt.plot(x[:plotlim], values1[0][:plotlim], label='Random start config')
    plt.plot(x[:plotlim], values2[0][:plotlim], label='Ordered start config')
    plt.xlabel("Monte Carlo cycles", fontsize=12)
    plt.ylabel(r'$\langle E\rangle$', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.grid()

    plt.subplot(222)
    plt.title(r'T = 2.4')
    plt.plot(x[:plotlim], values3[0][:plotlim], label='Random start config')
    plt.plot(x[:plotlim], values4[0][:plotlim], label='Ordered start config')
    plt.xlabel("Monte Carlo cycles", fontsize=12)
    plt.ylabel(r'$\langle E\rangle$', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.grid()

    plt.subplot(223)
    plt.title(r'T = 1.0')
    plt.plot(x[:plotlim], values1[1][:plotlim], label='Random start config')
    plt.plot(x[:plotlim], values2[1][:plotlim], label='Ordered start config')
    plt.xlabel("Monte Carlo cycles", fontsize=12)
    plt.ylabel(r'$\langle |M|\rangle$', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.grid()

    plt.subplot(224)
    plt.title(r'T = 2.4')
    plt.plot(x[:plotlim], values3[1][:plotlim], label='Random start config')
    plt.plot(x[:plotlim], values4[1][:plotlim], label='Ordered start config')
    plt.xlabel("Monte Carlo cycles", fontsize=12)
    plt.ylabel(r'$\langle |M|\rangle$', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.grid()
    plt.savefig('c_values.pdf')

    plt.figure(figsize=(8,11))
    plt.subplot(211)
    plt.plot(x[:plotlim], np.add.accumulate(values1[4][:plotlim]), label="T = 1. Random")
    plt.plot(x[:plotlim], np.add.accumulate(values2[4][:plotlim]), label= "T = 1. Ordered")
    plt.xlabel("Monte Carlo cycles", fontsize=14)
    plt.ylabel("Accepted configurations", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=14)
    plt.grid()

    plt.subplot(212)
    plt.plot(x[:plotlim], np.add.accumulate(values3[4][:plotlim]), label= "T = 2.4. Random")
    plt.plot(x[:plotlim], np.add.accumulate(values4[4][:plotlim]), label= "T = 2.4. Ordered")
    plt.xlabel("Monte Carlo cycles", fontsize=14)
    plt.ylabel("Accepted configurations", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=14)
    plt.grid()
    plt.savefig("c_acc.pdf")
    plt.show()

def ex_d():

    def gauss(x, fac, mu, sigma):
        ex = np.exp(-(x - mu)**2/(2*sigma**2))
        return fac*ex

    filenames = ["ex_d_1", "ex_d_2.4"]
    values = []
    probs = []

    for i in range(2):
        n, m, t = get_dim(filenames[i]+"_prob")
        values.append(read_file(filenames[i], 2, 7))
        probs.append(read_file(filenames[i]+"_prob", m, n))

    v1 = np.zeros(7)
    v2 = v1.copy()

    for i in range(7):
        v1[i] = values[0][i][0]
        v2[i] = values[1][i][0]

    sortindex1 = probs[0][0].argsort()
    sortindex2 = probs[1][0].argsort()

    x1 = probs[0][0][sortindex1]*20**2
    mu1, sigma1 = v1[0]*20**2, np.sqrt(v1[5]*20**2)
    x1 = np.linspace(x1[0] + (x1[0] - x1[-1]), x1[-1], 2*len(x1) + 1)

    y1 = np.zeros(2*len(probs[0][1]) + 1)
    y1[len(probs[0][1]):-1] = probs[0][1][sortindex1]
    y1[y1==0] = np.nan

    mu2, sigma2 = v2[0]*20**2, np.sqrt(v2[5]*20**2)
    x2 = probs[1][0][sortindex2]*20**2
    y2 = probs[1][1][sortindex2]

    popt2, pcov2 = curve_fit(gauss, x2, y2, p0=[1, mu2, sigma2])


    plt.figure(figsize=(7,12))
    plt.subplot(211)
    plt.title(r'T = 1.0', fontsize=15)
    plt.plot(x1, y1, '+', label=r'$P(E_i)$')
    plt.xlim(x1[0], x1[-1])
    plt.axvline(x = mu1, color='black', label=r'$\langle E\rangle$')
    plt.axvspan(mu1 - sigma1, mu1 + sigma1, alpha=0.2, color='red', label=r'$\langle E\rangle \pm \sigma_E$')
    plt.xlabel(r'E', fontsize=15)
    plt.ylabel(r'P(E)', fontsize=15)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.legend(fontsize=15)
    plt.grid(True)

    plt.subplot(212)
    plt.title(r'T = 2.4', fontsize=15)
    plt.plot(x2, y2, '+', label=r'$P(E_i)$')
    plt.xlabel(r'E', fontsize=15)
    plt.ylabel(r'P(E)', fontsize=15)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.axvline(x = mu2, color='black', label=r'$\langle E\rangle$')
    plt.axvspan(mu2 - sigma2, mu2 + sigma2, alpha=0.2, color='red', label=r'$\langle E\rangle \pm \sigma_E$')
    plt.plot(x2, gauss(x2, *popt2), 'r:', label='Gaussian fit')
    plt.legend(fontsize=13)
    plt.grid()
    plt.savefig('d.pdf')
    plt.show()

def ex_e():
    filename40 = "ex_e_40x40_0.01"
    filename60 = "ex_e_60x60_0.01"
    filename80 = "ex_e_80x80_0.01"
    filename100 = "ex_e_100x100_0.01"
    filename140 = "ex_e_140x140_0.01"
    values40 = read_file(filename40, 41, 7)
    values60 = read_file(filename60, 41, 7)
    values80 = read_file(filename80, 41, 7)
    values100 = read_file(filename100, 41, 7)
    values140 = read_file(filename140, 41, 7)

    values = [values40, values60, values80, values100, values140]
    labels = ["40x40", "60x60", "80x80", "100x100", "140x140"]
    colors = ["blue", "red", "green", "orange", "purple"]
    x = np.linspace(2.2, 2.6, 41)
    plotlim = 20

    plt.figure(figsize=(12,10))
    plt.subplot(221)
    for i, v in enumerate(values):
        plt.plot(x[:plotlim], v[0][:plotlim], label=labels[i], color=colors[i])
        plt.plot(x[:plotlim], v[0][:plotlim], 'o', markersize=3, color=colors[i])
    plt.ylabel(r'$\langle E\rangle$', fontsize=15)
    plt.xlabel('kT', fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.grid()

    plt.subplot(222)
    for i, v in enumerate(values):
        plt.plot(x[:plotlim], v[4][:plotlim], label=labels[i], color=colors[i])
        plt.plot(x[:plotlim], v[4][:plotlim], 'o', markersize=3, color=colors[i])
    plt.ylabel(r'$\langle |M|\rangle$', fontsize=15)
    plt.xlabel(r'kT', fontsize=15)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid()

    plt.subplot(223)
    for i, v in enumerate(values):
        plt.plot(x[:plotlim], v[5][:plotlim], label=labels[i], color=colors[i])
        plt.plot(x[:plotlim], v[5][:plotlim], 'o', markersize=3, color=colors[i])
    plt.axvline(x =  2/(np.log(1 + np.sqrt(2))), linestyle=':', color='red', label=r'$T_C \approx 2.269$')
    plt.ylabel(r'$C_V/Jk$', fontsize=15)
    plt.xlabel(r'kT', fontsize=15)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid()

    plt.subplot(224)
    for i, v in enumerate(values):
        plt.plot(x[:plotlim], v[6][:plotlim], label=labels[i], color=colors[i])
        plt.plot(x[:plotlim], v[6][:plotlim], 'o', markersize=3, color=colors[i])
    plt.axvline(x =  2/(np.log(1 + np.sqrt(2))), linestyle=':', color='red', label=r'$T_C \approx 2.269$')
    plt.ylabel(r'$\chi$', fontsize=15)
    plt.xlabel(r'kT', fontsize=15)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid()
    plt.savefig('f.pdf')
    plt.show()


def ex_f():
    filename40 = "ex_e_40x40_2.25-2.29"
    filename60 = "ex_e_60x60_2.25-2.29"
    filename80 = "ex_e_80x80_2.25-2.29"
    filename100 = "ex_e_100x100_2.25-2.29"
    filename140 = "ex_e_140x140_2.25-2.29"
    values40 = read_file(filename40, 28, 7)
    values60 = read_file(filename60, 28, 7)
    values80 = read_file(filename80, 28, 7)
    values100 = read_file(filename100, 28, 7)
    values140 = read_file(filename140, 28, 7)

    values = [values40, values60, values80, values100, values140]
    labels = ["40x40", "60x60", "80x80", "100x100", "140x140"]
    colors = ["blue", "red", "green", "orange", "purple"]
    L = [40, 60, 80, 100, 140]
    var = []
    tcL = []
    Tc = []


    x = np.linspace(2.25, 2.29, 28)
    window = 8
    winIn = window//2
    v = np.ones(int(window))/window

    for i in range(5):
        values[i][5] /= x**2
        values[i][6] /= x
        var.append(np.convolve(values[i][5], v, "same"))
        plt.plot(x[winIn:-winIn], var[i][winIn:-winIn], label=labels[i], color = colors[i])
        plt.plot(x[winIn:-winIn], values[i][5][winIn:-winIn], '+', color = colors[i], alpha=0.5)
        tcL.append(x[var[i].argmax()])
        plt.scatter(tcL[i], var[i].max(), color=colors[i])

    k = 0
    for i in range(5):
        for j in range(k, 5, 1):
            if j != i:
                Tc.append((tcL[i]*L[i] - tcL[j]*L[j])/(L[i] - L[j]))
        k += 1

    print(np.mean(Tc))
    print(np.std(Tc))

    plt.xlabel(r'kT', fontsize=12)
    plt.ylabel(r'$C_V/Jk$', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.grid()
    plt.savefig('e_CV.pdf')
    plt.show()

#ex_b()
#ex_c()
#ex_d()
#ex_e()
#ex_f()
