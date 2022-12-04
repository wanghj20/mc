import numpy as np
import matplotlib.pyplot as plt
import math as m
from scipy import integrate
from mpmath import quad
from scipy.stats import norm


def primRoots(modulo, generator=True):
    coprime_set = {num for num in range(1, modulo) if m.gcd(num, modulo) == 1}
    if generator:
        return (g for g in range(1, modulo) if coprime_set == {
            pow(g, powers, modulo) for powers in range(1, modulo)})

    return [g for g in range(1, modulo) if coprime_set == {
        pow(g, powers, modulo) for powers in range(1, modulo)}]

def Phi(x):
    return norm.cdf(x)

def PhiInv(x):
    return norm.ppf(x)

def inter1(x,y):
    return (x-y)*np.exp(2*a*PhiInv(x))*np.exp(PhiInv(x)**2/2)*np.sqrt(2*np.pi)

def inter2(x):
    return (x**2)*np.exp(2*a*PhiInv(x))*np.exp(PhiInv(x)**2/2)*np.sqrt(2*np.pi)

def theta1(y):
    theta_a1 = integrate.quad(lambda x: (x-y)*np.exp(2*a*PhiInv(x))*np.exp(PhiInv(x)**2/2)*np.sqrt(2*np.pi), y, 1/2)
    theta_a2 = integrate.quad(lambda x: (x**2)*np.exp(2*a*PhiInv(x))*np.exp(PhiInv(x)**2/2)*np.sqrt(2*np.pi), 0, 1/2)
    return 2 * (theta_a1[0]-theta_a2[0])
    #return theta_a1[1],theta_a2[1]

def theta(y):
    h = 2/m*np.log(np.pi*m)
    print(h)
    r1 = np.ones(m+1, dtype=object)
    print(inter1(v1(-3*h, 1/1009), 1/1009)*diffv1(-3*h,1/1009))

    for j in range(m+1):
        r1[j] = inter1(v1(-j*h,y),y) * diffv1(-j*h,y)

    print(r1[5])
    theta_a1 = h*np.sum(r1)
    r2 = np.ones(m+1, dtype=object)
    for j in range(1, m+1):
        r2[j] = inter2(v2(-j*h))*diffv2(-j*h)

    theta_a2 = h*np.sum(r2)
    return 2*(theta_a1 - theta_a2)


def v1(t, y):
    return (1/2-y)*np.tanh(np.pi/2*np.sinh(t)) + 1/2

def v2(t):
    return 1/2*np.tanh(np.pi/2*np.sinh(t)) + 1/2

def diffv1(t,y):
    return (1/2-y)*(1-np.tanh(np.pi/2*np.sinh(t))**2)*np.pi/2*np.cosh(t)

def diffv2(t):
    return 1/2*(1-np.tanh(np.pi/2*np.sinh(t))**2)*np.pi/2*np.cosh(t)

def fastpod(n,s,aa,omega,gamma,g):
    """Fast cbc algorithm with Pod weights."""
    z = np.ones(s)
    e2 = np.zeros(s)
    Gamma_r = np.ones(s)
    for j in range(s):
        Gamma_r[j] = (j+1) ** aa
    m1 = int((n-1)/2)
    E2 = np.zeros(m1, dtype=complex)
    perm = np.zeros(m1)
    perm[0] = 1
    for j in range(1, m1):
        perm[j] = (perm[j-1] * g) % n
    perm = np.minimum(perm, n-perm)
    print(perm)
    psi = np.zeros(m1)
    for t in range(m1):
        psi[t] = omega((perm[t])/n)
    psi0 = omega(0)
    fft_psi = np.fft.fft(psi)
    prev_e2 = 0
    q = np.zeros((m1, s))
    q[:, 0] = np.ones(m1)
    q0 = np.zeros(s)
    q0[0] = 1
    #E2 = np.real(np.fft.ifft(fft_psi * np.fft.fft(q@np.diag(Gamma_r).sum(axis=1))))
    #q[:, 0] = q[:, 0] + gamma[0]*Gamma_r[0] * np.r_[ psi[0:0:-1], psi[0], psi[m1:0:-1]] * q[:, -1]
    #prev_e2 = gamma[0] * ((psi0 * np.dot(q0, Gamma_r)) + 2 * np.min(E2)) / n
    #q0[0] = q0[0] + gamma[0] * Gamma_r[0] * psi0 * q0[-1]
    for k in range(s):
        alpha_r = Gamma_r[:k]
        diagalpha = np.diag(alpha_r)
        E2 = np.real(np.fft.ifft(fft_psi * np.fft.fft((q@np.diag(Gamma_r)).sum(axis=1))))
        print((q@np.diag(Gamma_r)).sum(axis=1))
        w = np.argmin(E2)
        z[k] = perm[w]
        e2[k] = prev_e2 + gamma[k] * ( psi0 * (np.dot(q0, Gamma_r)) + 2 * np.min(E2)) / n
        for l in range(k, 0, -1):
            q[:, l] = q[:, l] + gamma[k]*Gamma_r[l] * np.r_[ psi[w:0:-1], psi[0], psi[m1:w:-1]] * q[:, l-1]
        #print(q)
        for l in range(k, 0, -1):
            q0[l] = q0[l] + gamma[k] * Gamma_r[l] * psi0 * q0[l-1]
        prev_e2 = e2[k]

    return prev_e2, z

if __name__ == '__main__':
    yg = next(primRoots(1009))
    #print(yg)
    s = 7
    kappa = 0.01
    eta = 3.1
    lam = 0.51
    a = 1/4
    m = 1000
    gamma1 = np.ones(s)
    for i in range(s):
        gamma1[i] = (kappa/((i+1)**eta))**(1/(lam+1))
    #print(gamma1[0])
    #print(theta1(1/5),theta1(4/5))
    print(fastpod(n=1009, s=7, aa=2/(1+lam), omega=theta1, gamma=gamma1, g=yg))
    #print(theta1(20/1009))





