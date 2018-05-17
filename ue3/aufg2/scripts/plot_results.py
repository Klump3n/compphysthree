#!/usr/bin/env python3
"""
Plot the results...

"""
import numpy as np


def plot_vals(laplace_vals, vec_add_vals, vec_scale_vals, save_name=None):

    if save_name is not None:
        import matplotlib as mpl
        mpl.use('Agg')
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(laplace_vals, label='Laplace Produkt')
    plt.plot(vec_scale_vals, label='Vektor Skalierung')
    plt.plot(vec_add_vals, label='Vektor Addition')
    plt.yscale('log')
    plt.xlabel(r'Konfigurations ID/Gitterdimension')
    plt.ylabel(r'Speedup durch Verwendung der GPU')
    plt.title('Speedup von Operationen auf der GPU')
    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout()

    if save_name is None:
        plt.show()
    else:
        plt.savefig('../figures/{}.pdf'.format(save_name))

if __name__ == '__main__':
    factorizations = np.genfromtxt(
        'factorizations', delimiter=',',
        usecols=(0,1,2,3,4), dtype=int)
    laplace = np.genfromtxt(
        'laplace_speedup_results', delimiter=',',
        usecols=(0,1,2), dtype=float, skip_header=1)
    vec_add = np.genfromtxt(
        'vector_add_speedup_results', delimiter=',',
        usecols=(0,1,2), dtype=float, skip_header=1)
    vec_scale = np.genfromtxt(
        'vector_scale_speedup_results', delimiter=',',
        usecols=(0,1,2), dtype=float, skip_header=1)

    laplace_speedup = laplace[:, 2]
    vec_add_speedup = vec_add[:, 2]
    vec_scale_speedup = vec_scale[:, 2]

    plot_vals(laplace_speedup, vec_add_speedup, vec_scale_speedup, save_name='speedup')

    # for finding weird behaviour
    for it, val in enumerate(laplace_speedup):
        if it > 0:
            val_left = laplace_speedup[it-1]
            val = laplace_speedup[it]

            if (max(val, val_left) / min(val, val_left)) > 1.5:
                print('configuration = [N, gridDim.x, gridDim.y, blockDim.x, blockDim.y]')
                print('configuration {} has speedup {}'.format(factorizations[it-1], laplace_speedup[it-1]))
                print('configuration {} has speedup {}'.format(factorizations[it], laplace_speedup[it]))
                print('ratio is {}\n'.format((max(val, val_left) / min(val, val_left))))
