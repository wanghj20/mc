import numpy as np
import scipy.stats
from scipy.stats import qmc
from scipy import integrate

def f(x):
    ff = np.zeros(5)
    for i in range(5):
        ff[i] = A[i]*x[i]**2*(1 + 1/2*np.sin(sum(x) - x[i]))
    return np.exp(np.sum(ff))

def g(x):
    return np.exp(x**2)

def p(x):
    return np.exp(sum(A*x**2))/C

def cond(x,t):
	con = (N - 1)*(M - 1)*p(x)
	return con/(con + t*(M - p(x)))

def AR_RQMC(function,shift):
	sample = []
	i = 0 #number of samples rejected
	j = 0 #number of samples accepted
	k = 0 #total number of samples
	while j < N:
		x = np.modf(k*z/(M*N) + shift)[0]
		sample.append(x)
		x1 = x[0:5]
		x2 = x[-1]
		if x2 < function(x1)/M:
			j = j + 1
		else:
			i = i + 1
		k = k +1
	return sample, i

#constant
A = np.array([1,1/2,1/5,1/5,1/5])
#print(f([1,2,3,4,5]))
z = np.array([243.,  86.,  72., 165., 471,308])
m = 64
s = 6
N = 4001
C_1 = integrate.quad(g,0,np.sqrt(A[0]))[0]/np.sqrt(A[0])
C_2 = integrate.quad(g,0,np.sqrt(A[1]))[0]/np.sqrt(A[1])
C_3 = integrate.quad(g,0,np.sqrt(A[2]))[0]/np.sqrt(A[2])
C = C_1*C_2*C_3**3
M = 1/C*(np.exp(2.2))

#randomly shifted
np.random.seed(7)
delt = np.random.uniform(size=(m, s))
I_w = np.zeros(m)
for l in range(m):
    sam ,t = AR_RQMC(p,delt[l,:])
    samp = sam[0:N+t-1]
    samp1 = sam[-1]
    I_w[l] = (sum([cond(y[0:5],t)*f(y[0:5])/p(y[0:5]) for y in samp]) + f(samp1[0:5])/p(samp1[0:5])) / N


#standard deviation
es = np.mean(I_w)
sd = np.sqrt(sum((I_w-es)**2) / (m-1) / m)

print("conditional smooth：",[es, sd])


