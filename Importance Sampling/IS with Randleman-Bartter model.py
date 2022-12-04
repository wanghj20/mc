import math as m

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

from cbc import cbcsearch, B2, gammas


dim = 10
r0 = 0.1
sigma = 0.01
zs = -1*0.2*np.ones(dim)
N = 97
R = 15


def diff(z):
    r = np.zeros(dim)
    rs = np.zeros(dim)
    for i in range(dim):
        r[i] = r0 * np.exp(-i*sigma**2*(1/2) + sigma*np.sum(z[i:]))
    for i in range(dim):
        rs[i] = -1*(sigma**2*r[i])/(1 + r[i]**2)
    F = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            if i <= j:
                F[i, j] = np.sum(rs[j:])
            if i > j:
                F[i, j] = np.sum(rs[i:])
    return F


def FF(y):
    rk = np.zeros(dim)
    rsi = np.zeros(dim)
    ff = np.zeros(dim)
    for i in range(dim):
        rk[i] = r0 * np.exp(-i*sigma**2*(1/2) + sigma*np.sum(y[i:]))
    for i in range(dim):
        rsi[i] = sigma*rk[i]/(1+rk[i])
    for i in range(dim):
        ff[i] = -1*np.sum(rsi[i:])-y[i]
    return ff


def newton(x):
    x1 = np.copy(x)
    i = 0
    delta = np.copy(x)
    while np.linalg.norm(delta) > 1.e-5 and i < 50:
        x1 = x - np.linalg.inv(diff(x)-np.eye(dim))@FF(x)
        delta = x1 - x
        x = x1
        i = i + 1
    return x


def g(x):
    rk = np.zeros(dim)
    for i in range(dim):
        rk[i] = r0 * np.exp(-i*sigma**2*(1/2) + sigma*np.sum(x[i:]))
    gg = -1*np.sum(np.log(1+rk))
    return gg


def f(x):
    y = np.matmul(L_star, x) + z_star
    fx = m.pow(2*np.pi, dim/2)*m.exp(g(y) - 1/2*y.T@diff(z_star)@y - z_star.T@np.linalg.inv(sigma_star)@L_star@x - \
                                     1/2*z_star.T@np.linalg.inv(sigma_star)@z_star)
    return fx


def mserror(r, s, n):
    delt = np.random.uniform(size=(r, s))
    qq = np.zeros((n, r))
    for i in range(n):
        for j in range(r):
            qq[i, j] = f(norm.ppf(np.modf(i*z_optimal/n + delt[j, :])[0]))
    qmc_sn = np.zeros(r)
    for i in range(r):
        qmc_sn[i] = 1/n*np.sum(qq[:, i])
    qmc_rns = 1/r*np.sum(qmc_sn)
    mse = 1/(r*(r-1))*np.sum((qmc_rns - qmc_sn[i])**2)
    return mse


if __name__ == '__main__':
    gamma = gammas(dim)

    z_star = newton(zs)
    # print(z_star)
    sigma_star = -1*np.linalg.inv(diff(newton(zs)) - np.eye(dim))
    L_star = np.linalg.cholesky(sigma_star)

    z_optimal = cbcsearch(gamma, N, dim)

    print(mserror(R, dim, N))
