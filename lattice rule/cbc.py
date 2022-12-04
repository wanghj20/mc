import numpy as np
import matplotlib.pyplot as plt
import math as m

dim = 100
N = 4001

def gammas(d):
    ga = np.zeros(d)
    for i in range(d):
        ga[i] = 1/((i+1)**2)
    return ga


def B2(x):
    return x**2-x+1/6


def cbcsearch(gamma, n, d):
    
    omega = np.zeros((n-1, n))
    for z in range(n-1):
        for k in range(n):
            omega[z, k] = (2*np.pi**2)*B2(m.modf(k*(z+1)/n)[0])
    
    poly = np.ones(n)
    for i in range(n):
        poly[i] = 1 + gamma[0]*omega[0, i]
 
    z_opt = np.zeros((d,), dtype=int)
    z_opt[0] = 1
    for s in range(1, d):
        error2 = -1+1/n*(1+gamma[s]*omega)@poly
        #print(error2)
        z_opt[s] = np.argmin(error2)+1
        v = np.zeros(n-1)
        v[z_opt[s]-1] = 1
        poly = np.diagflat(1+gamma[s]*v@omega)@poly

  
    return z_opt


if __name__ == '__main__':
    gamma = gammas(dim)
    print(cbcsearch(gamma, N, dim))
