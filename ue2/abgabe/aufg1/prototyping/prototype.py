#!/usr/bin/env python3
"""
Prototype for that thaaaang.

"""
import numpy as np
from laplace import laplace


def bare_vec():
    return np.asarray(
        [
            10.30, 19.80, 10.50, 11.50, 8.10, 25.50, 7.40, 23.60,
            4.10, 20.50, 18.60, 17.10, 24.20, 25.10, 22.70, 7.00,
            12.40, 19.40, 8.40, 24.80, 2.70, 23.20, 23.10, 14.10,
            11.80, 9.00, 4.60, 9.90, 5.10, 15.90, 20.10, 15.40,
            10.20, 5.00, 1.30, 18.30, 4.90, 8.80, 16.30, 9.00,
            3.70, 9.30, 0.50, 2.30, 8.80, 23.30, 9.40, 21.20,
            17.10, 17.80, 20.50, 19.80, 15.50, 18.00, 8.40, 1.70,
            1.40, 13.00, 11.60, 6.50, 3.30, 6.10, 22.00, 13.50,
        ]
    )


def vec():
    bare_vector = bare_vec()
    reshaped_bare_vector = bare_vector.reshape((8, 8))
    embed = np.zeros((10, 10))
    embed[1:9,1:9] = reshaped_bare_vector
    return embed


def flatten_embed():
    return vec().flatten()


def manual():
    N = 8
    n = N+2
    npts = n*n

    b = np.zeros(n*n)

    x = flatten_embed()         # len 100 for N = 8

    for i in range(npts):
        for j in range(npts):

            # no boundary
            if (
                    (i % n == 0) or
                    (int(i / n) == 0) or
                    (i % n == n-1) or
                    (int(i / n) == n-1)
            ):
                continue

            # diagonal
            if (i == j):
                b[i] += -4 * x[j]

            # first off diagonal
            if (
                    ((i == j+1) and (i % n != 0))
                    or
                    ((i == j-1) and (j % n != 0))
            ):
                b[i] += 1 * x[j]

            # second off diagonal
            if ((i == j+n) or (i == j-n)):
                b[i] += 1 * x[j]

    return b


def extract_core(embedded_result):
    reshaped = embedded_result.reshape((10, 10))
    return reshaped[1:9, 1:9].flatten()


if __name__ == '__main__':

    # do the thing manually with boundary
    man = manual().reshape((10, 10))
    print(man)

    boundary_product = extract_core(man)

    # just multiiply
    bare_product = np.dot(-1*laplace(8), bare_vec())

    print(bare_product.reshape((8, 8)))

    # compare results
    print('A*x == the manual thing? {}'.format(np.allclose(boundary_product, bare_product)))
