import numpy as np

En = -7.98393
E2n = 63.8714
Mn = 3.99454
M2n = 15.9732
CVn = 0.128329
Sn = 0.016043

Ea = -7.9837
E2a = 63.8696
Ma = 3.99456
M2a = 15.9728
CVa = 0.13016
Sa = 0.0162916

relE = abs(En - Ea)/Ea
relE2 = abs(E2n - E2a)/E2a
relM = abs(Mn - Ma)/Ma
relM2 = abs(M2n - M2a)/M2a
relCV = abs(CVn - CVa)/CVa
relS = abs(Sn - Sa)/Sa


rel = [relE, relE2, relM, relM2, relCV, relS]

for i in rel:
    print(i)
