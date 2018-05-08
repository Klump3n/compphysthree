#!/usr/bin/env python3
"""
Implementation of the conjugate gradient method.

"""
import time
import numpy as np
from laplace import laplace

from cg import cg, two_norm
from scipy.sparse.linalg import cg as scicg  # for comparison


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


def plot_residuum(r, title=None, save_name=None):
    """
    Plot for the errors.

    """
    if save_name is not None:
        import matplotlib as mpl
        mpl.use('Agg')
    import matplotlib.pyplot as plt

    # normalize to r_0
    r_normalised = r / r[0]

    plt.figure()
    plt.plot(r_normalised)
    plt.yscale('log')
    plt.ylabel(r'Relatives Residuum $\| r_{k} \| / \| r_{0} \|$')
    plt.xlabel(r'Iterationsschritt $k$')
    if title:
        plt.title('Residuum: {}'.format(title))
    plt.grid()
    plt.tight_layout()

    if save_name is None:
        plt.show()
    else:
        plt.savefig('../figures/residuum_{}.pdf'.format(save_name))


def plot_errors(e, title=None, save_name=None):
    """
    Plot for the errors.

    """
    if save_name is not None:
        import matplotlib as mpl
        mpl.use('Agg')
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(e)
    plt.yscale('log')
    # plt.ylabel(r'Fehler $\| r_{k} \| / \| r_{0} \|$')
    plt.xlabel(r'Iterationsschritt $k$')
    if title:
        plt.title('Fehler: {}'.format(title))
    plt.grid()
    plt.tight_layout()

    if save_name is None:
        plt.show()
    else:
        plt.savefig('../figures/errors_{}.pdf'.format(save_name))


def plot_residuum_and_errors(r, e, title=None, save_name=None):
    """
    Plot for the errors.

    """
    if save_name is not None:
        import matplotlib as mpl
        mpl.use('Agg')
    import matplotlib.pyplot as plt

    # normalize to r_0
    r_normalised = r / r[0]

    plt.figure()
    plt.plot(r_normalised, label='Residuum')
    plt.plot(e, label='Fehler')
    plt.yscale('log')
    # plt.ylabel(r'Relatives Residuum $\| r_{k} \| / \| r_{0} \|$')
    plt.xlabel(r'Iterationsschritt $k$')
    if title:
        plt.title('Residuum und Fehler: {}'.format(title))
    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout()

    if save_name is None:
        plt.show()
    else:
        plt.savefig('../figures/residuum_and_errors_{}.pdf'.format(save_name))


if __name__ == '__main__':
    """
    Standalone call.

    """
    N = 64                            # lattice size
    A = -1 * laplace(int(np.sqrt(N)))  # solve the discretized laplace eqn

    print('Problem size is {}x{}\n'.format(N, N))

    A_eigvals, A_eigvecs = np.linalg.eig(A)
    A_sorted_indices = np.abs(A_eigvals).argsort()
    A_smallest_eigvec = A_eigvecs[:, A_sorted_indices[0]]  # smallest eigenvector will be first

    unit_vector = np.zeros(N)
    unit_vector[np.random.randint(N)] = 1  # some random entry is set to one

    ones_vector = np.dot(A, np.ones(N))

    sum_of_eigvecs = A_eigvecs[:, np.random.randint(N)] + A_eigvecs[:, np.random.randint(N)]

    b_dict = {
        'random': {
            'b_value': np.random.rand(N),
            'title': r'Zufälliger Vektor',
            'save_name': 'zufaelliger_vektor'
        },
        'smallest_eigenvec': {
            'b_value': A_smallest_eigvec,
            'title': r'Kleinster Eigenvektor',
            'save_name': 'kleinster_eigenvektor'
        },
        'unit_vector': {
            'b_value': unit_vector,
            'title': r'Einheitsvektor',
            'save_name': 'einheitsvektor'
        },
        'ones': {
            'b_value': ones_vector,
            'title': r'Einsvektor',
            'save_name': 'einsvektor'
        },
        'sum': {
            'b_value': sum_of_eigvecs,
            'title': r'Summe zweier Eigenvektoren',
            'save_name': 'summe_eigenvektoren'
        }
    }

    # for every vector ...
    for pair in b_dict:

        # ... extract values from dict ...
        b_value = b_dict[pair]['b_value']
        plot_title = b_dict[pair]['title']
        save_name = b_dict[pair]['save_name']

        print(plot_title)

        # ... and do some magic
        t1 = time.monotonic()
        x_0, R = cg(A, b_value)
        t2 = time.monotonic()
        print('{:.4f} seconds with own cg'.format(t2 - t1))

        # get the lengths of the residuum and the error
        R_norm = calc_residuum(R)
        E = calc_errors(R, A)

        # plot the results
        plot_residuum(R_norm, title=plot_title, save_name=save_name)
        plot_errors(E, title=plot_title, save_name=save_name)
        plot_residuum_and_errors(R_norm, E, title=plot_title, save_name=save_name)

        # # for comparison
        # t1 = time.monotonic()
        # x_1 = np.linalg.solve(A, b_value)  # What I DONT understand is why its sometimes fast and sometimes slow...
        # t2 = time.monotonic()
        # print('{:.4f} seconds with linalg solver'.format(t2 - t1))

        # # for comparison
        # t1 = time.monotonic()
        # x_2 = scicg(A, b_value)
        # t2 = time.monotonic()
        # print('{:.4f} seconds with scipy conjugate gradient solver'.format(t2 - t1))

        # courtesy newline
        print()
