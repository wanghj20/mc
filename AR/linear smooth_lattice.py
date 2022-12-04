import numpy as np
from scipy.stats import norm
from scipy import integrate

def f(x):
    ff = np.zeros(5)
    for i in range(5):
        ff[i] = A[i]*x[i]**2*(1+1/2*np.sin(sum(x)-x[i]))
    return np.exp(np.sum(ff))
def g(x):
    return np.exp(x**2)
def p(x):
    return np.exp(sum(A*x**2))/C
def a(x):
    return np.exp(sum(A*x**2))/C - 1/2*L
def b(x):
    return np.exp(sum(A*x**2))/C + 1/2*L

def Wi_linear(x):
    if L*x[5]<=a(x[0:5]):
        return 1
    elif a(x[0:5])<=L*x[5]<=p(x[0:5]):
        W = 1+(p(x[0:5])-b(x[0:5]))/(p(x[0:5])-a(x[0:5]))/(b(x[0:5])-a(x[0:5]))*(L*x[5]-a(x[0:5]))
        return W
    elif p(x[0:5])<=L*x[5]<=b(x[0:5]):
        W = (p(x[0:5])-a(x[0:5]))/(p(x[0:5])-b(x[0:5]))/(b(x[0:5])-a(x[0:5]))*(L*x[5]-b(x[0:5]))
        return W
    elif b(x[0:5])<=L*x[5]<=L:
        return 0

def weight_sample_SLR(function,shift,N):
    sample = []
    sum = 0
    i = 0
    while sum < N:
        x = np.modf(i*z/N + shift)[0]
        sample.append(x)
        sum = sum + function(x)
        i = i+1
    return sample


A = np.array([1,1/2,1/5,1/5,1/5])
#print(f([1,2,3,4,5]))
z = np.array([243.,  86.,  72., 165., 471,308])
m = 20
s = 6
N = 4001
C_1 = integrate.quad(g,0,np.sqrt(A[0]))[0]/np.sqrt(A[0])
C_2 = integrate.quad(g,0,np.sqrt(A[1]))[0]/np.sqrt(A[1])
C_3 = integrate.quad(g,0,np.sqrt(A[2]))[0]/np.sqrt(A[2])
C = C_1*C_2*C_3**3
L = 1/C*(np.exp(3))

#randomly shifted
np.random.seed(3)
delt = np.random.uniform(size=(m, s))
I_w = np.zeros(m)
for i in range(m):
    sample = weight_sample_SLR(Wi_linear,delt[i,:],N)
    I_w[i] = sum([Wi_linear(y)*f(y[0:5])/p(y[0:5]) for y in sample]) / N

es = np.mean(I_w)
var = sum((I_w-es)**2) / (m-1) / m

print("linear：",[es, var])



