# -*- coding: utf-8 -*-
import numpy as np


def scramble(A, p=16, base=2):
    """Nested uniform scramble in base 2.

    A: (n, d) matrix with entries in [0, 1]
    p: precision parameter

    A is the matrix consists of n d-dimensional points. Each row of A is a
    point and each column is a dimension. The permutations used to scramble
    are identical within one dimension, and independent between different
    dimensions.
    """
    n, d = A.shape
    X = np.zeros((n, d))
    M = 1 << p

    for j in range(d):
        size = 1
        permutation = np.random.randint(2, size=size)
        for _ in range(p-1):
            size *= 2
            permutation = permutation.repeat(2)*2 + np.random.randint(2,
                size=size)

        def g(a):
            k = int(a) # integer part
            if k == M:
                k -= 1
            f = a - k # fractional part

            pi = permutation[k>>1]

            return ((k^pi)+f) / M

        X[:, j] = list(map(g, A[:, j]*M))

    return X


if __name__ == '__main__':
    pass
    d = 6
    m = 2
    p = 17
    seed = 110

    import scipy.stats

    np.random.seed(seed)
    A = scipy.stats.qmc.Sobol(d=d, scramble=False).random(1)

    X = scramble(A, p)
    print(X)

