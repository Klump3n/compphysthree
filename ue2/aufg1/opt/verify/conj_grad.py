#!/usr/bin/env python3
"""
Implementation of the conjugate gradient method.

"""
import numpy as np
from laplace import laplace

from cg import cg
from prototype import bare_vec


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
