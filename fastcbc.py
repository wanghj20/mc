import numpy as np
import matplotlib.pyplot as plt
import math as m
import scipy as sci
from scipy.stats import norm


# def pai(g, r):
#     z = np.zeros((1,N))
#     z[0] = g[0]
#     for i in range(1, N):
#         k = r**(i-1) % N
#         z[i] = g[k]
#     return z

def B2(x):
    return x**2-x+1/6

#求本原根算法

def primRoots(modulo, generator=True):
    """Find all primitive roots for the given modulo."""
    coprime_set = {num for num in range(1, modulo) if m.gcd(num, modulo) == 1}

    if generator:
        return (g for g in range(1, modulo) if coprime_set == {
            pow(g, powers, modulo) for powers in range(1, modulo)})

    return [g for g in range(1, modulo) if coprime_set == {
        pow(g, powers, modulo) for powers in range(1, modulo)}]

#循环矩阵的写法
# def circle(w,n):
#     circle1 = np.diagflat(w)
#     for i in range(1, n-1):
#         circle1 = circle1 + np.diagflat(w[n-1-i]*np.ones(n-1-i),i)
#     for i in range(1, n-1):
#         circle1 = circle1 + np.diagflat(w[i]*np.ones(n-1+i),-i)
#     return circle1



def fastrank1(n,s_max,gamma,g,omega=B2):
    z_opt = np.zeros(s_max)
    n1 = int((n-1)/2)
    perm = np.zeros(n1, 'i')
    perm[0] = 1
    perm = np.minimum(perm, n-perm)
    for i in range(1, n1):
        perm[i] = (perm[i-1]*g)%n
    psi0 = omega(0)
    psi = omega(perm/n)
    fft_psi = np.fft.fft(psi)
    q = np.ones(n1)
    q0 = 1
    for i in range(1, s_max):
        E2 = np.fft.ifft(np.diagflat(fft_psi) @ np.fft.fft(q))
        w = np.argmin(E2)
        z_opt[i] = perm[w]
        q = (1 + gamma[i]) * np.diagflat(psi) * q
        q0 = (1 + gamma[i] * psi0) * q0
    min_E2 = np.min(E2)
    e2 = -1 + 1/n*(q0 + 2*np.sum(q) + gamma[s_max]*(psi0*q0 + 2*min_E2))
    return e2, z_opt

if __name__ == '__main__':
    N = 4001
    yg = next(primRoots(N))
    s = 100
    gammas = 1/s*np.ones(s)
    print(fastrank1(N, s, gammas, yg, omega=B2))
    print()












# def fastsearch(gamma, n, d):
#     omega = np.zeros((n-1, n))
#     for z in range(n-1):
#         for k in range(n):
#             omega[z, k] = (2*np.pi**2)*B2(m.modf(k*(z+1)/n)[0])
#     poly = np.ones(n)
#     for i in range(n):
#         poly[i] = 1 + gamma[0]*omega[0, i]
#     q = pai(poly, x)
#     q1 = pai(poly, x)[1:n]
#     b = omega[0, 0]*np.ones(n)
#     G = np.concatenate((b, omega), axis=0)
#     phi = pai(G[:, 1], x)
#     pair = G[:, 1]@np.linalg.inv(phi)
#     phi1 = np.fft.fft(pai(G[:, 1], x)[1:n])
#     w_opt = np.zeros((d,), dtype=int)
#     w_opt[0] = 1
#     for s in range(1, d):
#         Error2 = 1/n * (gamma[s] * (np.fft.ifft(np.diagflat(phi1 @ np.fft(q1))) + gamma[s] * phi[0] * q[0] + np.sum(q)))-1
#         w_opt[s] = np.argmin(Error2)+1
#         v = np.zeros(n-1)
#         v[w_opt[s]-1] = 1
#         poly = np.diagflat(1+gamma[s]*v@np.concatenate((phi[0]*np.ones(n-1).T,circle(phi1,n)),axis=1))@poly
#     z_opt = pair@w_opt
#     return z_opt






# a = np.zeros(6)
# for i in range(0, 5):
#     a[i] = i
#
# print(1+a)
# A = np.linspace(0.001,0.999,100)
#
# plt.plot(A, norm.ppf(A))
# plt.show()
# x2 = lambda x: m.exp(-1*(x - miu)**2/(2*si**2))*1/(m.sqrt(2*m.pi*si))
# x0 = sci.integrate.quad(x2, 1, np.inf)[0]
# print(x0)
# #print(sci.integrate.quad(x2, index[1], np.inf))
# y = np.zeros(10)
# for i in range(10):
#     y[i] = 1 - sci.integrate.quad(x2, index[i], np.inf)[0]
#
# xnew = np.random.uniform(0,1,10)
# for i in range(10):
#     z[i] = interpolate.interp1d(y, index)(xnew[i])
# print(z)
#
# f = interpolate.interp1d(z, xnew)
# print(f(0.4))
# plt.plot(index, y, 'o', xnew, ynew, '-')
# plt.show()

