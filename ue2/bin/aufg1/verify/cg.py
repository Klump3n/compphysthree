#!/usr/bin/env python3
"""
CG method.

"""
import numpy as np


def two_norm(x):
    """
    The two norm.

    """
    return np.linalg.norm(x, 2)


def cg(A, b):
    """
    CG function.

    """
    n = len(b)
    eps = np.finfo(float).eps   # machine accuracy
    tol = n * eps               # convergence criterion

    x = np.zeros_like(b)
    r = b - np.dot(A, x)
    p = r

    r_norm_0 = two_norm(r)

    R = []
    R.append(r)

    rr = np.dot(r, r)

    for k in range(n):          # we expect converged results after n steps
        Ap = np.dot(A, p)
        alpha = rr / np.dot(p, Ap)
        x = x + alpha * p
        r = r - alpha * Ap

        R.append(r)

        crit = two_norm(r) / r_norm_0
        if crit < tol:
            print(
                'Tolerance reached, have || r || / || r_0 || = {} for k = {}'.format(
                    two_norm(r), k+1))
            break

        rr_new = np.dot(r, r)
        p = r + rr_new / rr * p
        rr = np.copy(rr_new)

    else:
        print('Did not converge quickly enough after {} steps'.format(k+1))

    return x, R
