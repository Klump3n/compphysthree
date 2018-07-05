#!/usr/bin/env python3
"""
Get metrics.

"""
import os
from subprocess import check_output
import matplotlib.pyplot as plt
import scipy.optimize as so
import re
import numpy as np
import time

cwd = os.getcwd()

h = 0.0
lamb = 1.0
L = 4

mean_k = .2333
delta_k = .1
kappa = .2333

res = {}

def func(kappa, factor, kappa_c):
    return factor * np.sqrt(kappa - kappa_c)


def ploti(thing):

    plt.figure()
    plt.xlabel(r'$\kappa$')
    plt.ylabel(r'$|M|^2$')

    for index, L in enumerate([4, 8, 16]):
        Lvals = thing[L]

        x = []
        y = []
        y_err = []

        for kappa in Lvals:
            y.append(Lvals[kappa]['mean'])
            y_err.append(Lvals[kappa]['std'])
            x.append(kappa)

        kappa = []
        mean = []
        std = []

        for index, value in enumerate(y):

            if (value > .2):
                kappa.append(x[index])
                mean.append(y[index])
                std.append(y_err[index])

        popt, pcov = so.curve_fit(func, kappa, mean, p0=[4., .2], sigma=std)

        kappa_c = popt[1]
        kappa_end = kappa[-1] + .1
        kappa_linspace = np.linspace(kappa_c, kappa_end)

        print('{}'.format(L))
        print(r'pref = {} $\pm$ {}'.format(popt[0], np.sqrt(pcov[0,0])))
        print(r'kappa_c = {} $\pm$ {}'.format(popt[1], np.sqrt(pcov[1,1])))

        ebarlines, _, _ = plt.errorbar(x, y, yerr=y_err, label='L = {}'.format(L), capsize=2, ls='None', marker='_')
        pcolor = ebarlines.get_color()
        plt.plot(kappa_linspace, func(kappa_linspace, *popt), label=r'L = {}, '.format(L) + r'$\kappa_{c} = (' + '{:.4f}'.format(popt[1]) + r' \pm {:.4f})$'.format(np.sqrt(pcov[1,1])), color=pcolor, alpha=.65)

    plt.legend(loc='best')
    plt.grid()
    plt.tight_layout()
    plt.savefig('mag_over_kappa.pdf')


for L in [4, 8, 16]:
    res[L] = {}

    for k in range(-10, 11):
        kappa = mean_k * (1 + delta_k * k)

        result = check_output(['{}/../bin/mc_mag'.format(cwd) ,'{}'.format(lamb), '{}'.format(kappa), '{}'.format(h), '{}'.format(0.0), '{}'.format(2000), '{}'.format(L), '{}'.format(L), '{}'.format(L)])
        val = (re.search(r'1000\\t.*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+)', str(result)).groups(0))

        mean = float(val[1])
        std = float(val[2])

        print(kappa, mean, std)

        res[L][kappa] = {}
        res[L][kappa]['mean'] = mean
        res[L][kappa]['std'] = std

        time.sleep(1)

ploti(res)










# configs = {
#     '1': {'N': 2048, 'dimx': 128, 'dimy': 1},
#          }

# results = []

# for conf in config:

#     N = config[conf]['N']
#     dimx = config[conf]['dimx']
#     dimy = config[conf]['dimy']

#     speedup_array = []
#     cpu_time_array = []
#     gpu_time_array = []

#     for _ in range(10):
#         result = check_output(['{}/../bin/run-cg'.format(cwd) ,'{}'.format(N), '{}'.format(dimx), '{}'.format(dimy)])

#         speedup = (re.search('Speedup = ([0-9]+\.[0-9]+)', str(result)).groups(0)[0])
#         cpu_time = (re.search('cpu cg      elapsed ([0-9]+\.[0-9]+)', str(result)).groups(0)[0])
#         gpu_time = (re.search('gpu cg      elapsed ([0-9]+\.[0-9]+)', str(result)).groups(0)[0])

#         speedup_array.append(float(speedup))
#         cpu_time_array.append(float(cpu_time))
#         gpu_time_array.append(float(gpu_time))


#     print(N, dimx, dimy, np.mean(speedup_array))
#     results.append([N, dimx, dimy, np.mean(speedup_array)])


# print('')
# with open('{}/results_file'.format(cwd), 'w') as f:
#     f.write('N, blockX, blockY, avg_time\n')
#     results = sorted(results)
#     for res in results:
#         f.write('{}, \t{}, \t{}, \t{}\n'.format(res[0], res[1], res[2], res[3]))
