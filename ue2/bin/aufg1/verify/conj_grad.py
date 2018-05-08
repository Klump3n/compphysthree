#!/usr/bin/env python3
"""
Implementation of the conjugate gradient method.

"""
import time
import numpy as np
from laplace import laplace

from cg import cg, two_norm

from prototype import bare_vec

def calc_residuum(R):
    """
    For every vector in R calculate the normalized residuum r.

    """
    return_r = []
    for r in R:
        return_r.append(two_norm(r))

    return return_r

def calc_errors(R, A):
    """
    Calculate the errors for every r in R.

    """
    Ainv = np.linalg.inv(A)
    return_e = []
    for r in R:
        return_e.append(two_norm(np.dot(Ainv, r)))

    return return_e


if __name__ == '__main__':
    """
    Standalone call.

    """
    N = 64                            # lattice size
    A = -1 * laplace(int(np.sqrt(N)))  # solve the discretized laplace eqn

    b_value = bare_vec()
    x_0, R = cg(A, b_value)

    np.set_printoptions(precision=2)  # only show two digits
    print(x_0.reshape((8, 8)))
