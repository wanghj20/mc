import math as m
import numpy as np
from scipy.stats import qmc
from scipy import integrate

def f(x):
    return np.exp(x*A)
def p(x):
    return np.exp(sum(A*x**2))/C
def a(x):
    return (1+sum(A*x**2))/C
def b(x):
    return np.exp(sum(A))/C

def Wi_cubic(x):
    if L*x[-1]<=a(x[0:5]):
        return 1
    elif a(x[0:5])<=L*x[-1]<=p(x[0:5]):
        W = (b(x[0:5])-p(x[0:5]))/(2*(b(x[0:5])-a(x[0:5])))*(((2*L*x[-1]-a(x[0:5])-p(x[0:5]))/(p(x[0:5])-a(x[0:5])))**3-3*((2*L*x[-1]-a(x[0:5])-p(x[0:5])) \
            /(p(x[0:5])-a(x[0:5]))))+(p(x[0:5])+b(x[0:5])-2*a(x[0:5]))/(2*(b(x[0:5])-a(x[0:5])))
        return W
    elif p(x[0:5])<=L*x[-1]<=b(x[0:5]):
        W = (p(x[0:5])-a(x[0:5]))/(2*(b(x[0:5])-a(x[0:5])))*(((2*L*x[-1]-b(x[0:5])-p(x[0:5]))/(b(x[0:5])-p(x[0:5])))**3-3*((2*L*x[-1]-b(x[0:5])-p(x[0:5])) \
            /(b(x[0:5])-p(x[0:5]))))+(p(x[0:5])-a(x[0:5]))/(2*(b(x[0:5])-a(x[0:5])))
        return W
    elif b(x[0:5])<=L*x[-1]<=L:
        retur

def weight_sample_RQMC(function,N):
    sample = []
    sampler = qmc.Sobol(d, scramble=True)
    sum = 0
    while sum < N:
        x = sampler.random(1)[0]
        sample.append(x)
        sum = sum + function(x)
    return sample

#RQMC

M = 20
d = 6
N = 1024
C = 3
A = np.array([1,1/2,1/5,1/5,1/5])
C_1 = integrate.quad(f,0,np.sqrt(A[0]))[0]/np.sqrt(A[0])
C_2 = integrate.quad(f,0,np.sqrt(A[1]))[0]/np.sqrt(A[1])
C_3 = integrate.quad(f,0,np.sqrt(A[2]))[0]/np.sqrt(A[2])
C = C_1*C_2*C_3**3
L = 2/C*(m.exp(3))



I_w = np.zeros(M)
for i in range(M):
    sample = weight_sample_RQMC(Wi_cubic,N)
    I_w[i] = sum([Wi_cubic(y)*f(y[0:5])/p(y[0:5]) for y in sample]) / N



es = np.mean(I_w)
error = np.abs(es-1)
var = sum((I_w-es)**2) / (M-1) / M

print("cubic curve smooth：",[es, error, var])
