#!/usr/bin/env python3
"""
Get metrics.

"""
import os
from subprocess import check_output
import re
import numpy as np


N = 128
dimx = 256
dimy = 1

cwd = os.getcwd()

h = 0.0
lamb = 1.0
L = 4

mean_k = .2333
delta_k = .02
kappa = .2333

result = check_output(['{}/../bin/mc_mag'.format(cwd) ,'{}'.format(lamb), '{}'.format(kappa), '{}'.format(h), '{}'.format(0.0), '{}'.format(1000), '{}'.format(L), '{}'.format(L), '{}'.format(L)])

print(result)
# for L in [4, 8, 16]:
#     for k in range(-10, 11):
#         # kappa = mean_k + delta_k * k
#         kappa = mean_k * (1 + delta_k * k)


        # WAIT 1 SEC













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
