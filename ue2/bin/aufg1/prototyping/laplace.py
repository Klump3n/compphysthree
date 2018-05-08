#!/usr/bin/env python3
"""
Port of laplace.m by B. Leder.

"""
import numpy as np
import scipy.sparse


def laplace(N):
    """
    Diskreter Laplace-Operator (\Delta) in 2 Dimensionen

    e = ones(N,1);
    E = spdiags(e,0,N,N);
    F = spdiags([-1*e 2*e -1*e],[-1 0 1],N,N);
    A = kron(F,E)+kron(E,F);

    A=full(A);

    """
    e = np.ones(N)
    diags = np.asarray([0])
    E = scipy.sparse.spdiags(e, diags, N, N)
    F = scipy.sparse.spdiags(np.array([-1*e, 2*e, -1*e]), np.array([-1, 0, 1]), N, N)
    A = scipy.sparse.kron(F, E) + scipy.sparse.kron(E, F)
    return A.toarray()


if __name__ == '__main__':
    print(laplace(3))
