#!/usr/bin/env python3
"""
Get metrics.

"""
import os
from subprocess import check_output
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

def ploti(thing):
    import matplotlib.pyplot as plt

    plt.figure()
    plt.xlabel(r'$\kappa$')
    plt.ylabel(r'$|M|^2$')

    for L in thing:
        Lvals = thing[L]

        x = []
        y = []
        y_err = []

        for kappa in Lvals:
            y.append(Lvals[kappa]['mean'])
            y_err.append(Lvals[kappa]['std'])
            x.append(float(kappa))

        plt.errorbar(x, y, yerr=y_err, label='{}'.format(L), capsize=1, fmt='+')

    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()


for L in [4, 8, 16]:
    res[str(L)] = {}

    for k in range(-10, 11):
        kappa = mean_k * (1 + delta_k * k)

        result = check_output(['{}/../bin/mc_mag'.format(cwd) ,'{}'.format(lamb), '{}'.format(kappa), '{}'.format(h), '{}'.format(0.0), '{}'.format(2000), '{}'.format(L), '{}'.format(L), '{}'.format(L)])
        val = (re.search(r'1000\\t.*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+).*([0-9]+\.[0-9]+)', str(result)).groups(0))

        mean = float(val[1])
        std = float(val[2])

        print(kappa, mean, std)

        res[str(L)][str(kappa)] = {}
        res[str(L)][str(kappa)]['mean'] = mean
        res[str(L)][str(kappa)]['std'] = std

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
