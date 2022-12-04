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

def AR_RQMC(function):
	sample = []
	i = 0 #拒绝样本的点数
	j = 0 #接受样本的点数
	while j < N:
		x = np.random.uniform(size=s)
		sample.append(x)
		x1 = x[0:5]
		x2 = x[-1]
		if x2 < function(x1)/M:
			j = j + 1
		else:
			i = i + 1
	return sample, i

#常数设定
A = np.array([1,1/2,1/5,1/5,1/5])
#print(f([1,2,3,4,5]))
m = 64
s = 6
N = 4096
C_1 = integrate.quad(g,0,np.sqrt(A[0]))[0]/np.sqrt(A[0])
C_2 = integrate.quad(g,0,np.sqrt(A[1]))[0]/np.sqrt(A[1])
C_3 = integrate.quad(g,0,np.sqrt(A[2]))[0]/np.sqrt(A[2])
C = C_1*C_2*C_3**3
M = 1/C*(np.exp(2.2))


np.random.seed(7)
I_w = np.zeros(m)
for l in range(m):
    sam ,t = AR_RQMC(p)
    samp = sam[0:N+t-1]
    samp1 = sam[-1]
    I_w[l] = (sum([cond(y[0:5],t)*f(y[0:5])/p(y[0:5]) for y in samp]) + f(samp1[0:5])/p(samp1[0:5])) / N


#误差计算
#误差计算
es = np.mean(I_w)
sd = np.sqrt(sum((I_w-es)**2) / (m-1) / m)

print("条件光滑：",[es, sd])



