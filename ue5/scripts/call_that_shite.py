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

old_config = {
    '32_1': {'N': 32, 'dimx': 128, 'dimy': 1},
    '32_2': {'N': 32, 'dimx': 256, 'dimy': 1},

    '64_1': {'N': 64, 'dimx': 256, 'dimy': 1},
    '64_2': {'N': 64, 'dimx': 256, 'dimy': 2},
    '64_3': {'N': 64, 'dimx': 256, 'dimy': 4},
    '64_4': {'N': 64, 'dimx': 128, 'dimy': 2},
    '64_5': {'N': 64, 'dimx': 128, 'dimy': 4},
    '64_6': {'N': 64, 'dimx': 64, 'dimy': 4},
    '64_7': {'N': 64, 'dimx': 32, 'dimy': 4},
    '64_8': {'N': 64, 'dimx': 32, 'dimy': 1},

    '128_1': {'N': 128, 'dimx': 256, 'dimy': 1},
    '128_2': {'N': 128, 'dimx': 256, 'dimy': 2},
    '128_3': {'N': 128, 'dimx': 256, 'dimy': 4},
    '128_4': {'N': 128, 'dimx': 128, 'dimy': 4},
    '128_5': {'N': 128, 'dimx': 64, 'dimy': 4},
    '128_6': {'N': 128, 'dimx': 32, 'dimy': 8},
    '128_7': {'N': 128, 'dimx': 32, 'dimy': 4},

    '256_1': {'N': 256, 'dimx': 8, 'dimy': 8},
    '256_2': {'N': 256, 'dimx': 16, 'dimy': 8},
    '256_3': {'N': 256, 'dimx': 32, 'dimy': 4},
    '256_4': {'N': 256, 'dimx': 64, 'dimy': 4},
    '256_5': {'N': 256, 'dimx': 128, 'dimy': 4},
    '256_6': {'N': 256, 'dimx': 256, 'dimy': 4},

    '512_1': {'N': 512, 'dimx': 256, 'dimy': 1},
    '512_2': {'N': 512, 'dimx': 256, 'dimy': 4},
    '512_3': {'N': 512, 'dimx': 128, 'dimy': 4},
    '512_4': {'N': 512, 'dimx': 64, 'dimy': 8},

    '1024_1': {'N': 1024, 'dimx': 256, 'dimy': 1},
    '1024_2': {'N': 1024, 'dimx': 256, 'dimy': 2},
    '1024_3': {'N': 1024, 'dimx': 256, 'dimy': 4},
    '1024_4': {'N': 1024, 'dimx': 128, 'dimy': 8},
    '1024_5': {'N': 1024, 'dimx': 128, 'dimy': 4},
    '1024_6': {'N': 1024, 'dimx': 64, 'dimy': 4},
    '1024_7': {'N': 1024, 'dimx': 32, 'dimy': 4},

}

config = { 
    '2048_1': {'N': 2048, 'dimx': 128, 'dimy': 1},
         }

results = []

for conf in config:

    N = config[conf]['N']
    dimx = config[conf]['dimx']
    dimy = config[conf]['dimy']

    speedup_array = []
    cpu_time_array = []
    gpu_time_array = []

    for _ in range(10):
        result = check_output(['{}/../bin/run-cg'.format(cwd) ,'{}'.format(N), '{}'.format(dimx), '{}'.format(dimy)])

        speedup = (re.search('Speedup = ([0-9]+\.[0-9]+)', str(result)).groups(0)[0])
        cpu_time = (re.search('cpu cg      elapsed ([0-9]+\.[0-9]+)', str(result)).groups(0)[0])
        gpu_time = (re.search('gpu cg      elapsed ([0-9]+\.[0-9]+)', str(result)).groups(0)[0])

        speedup_array.append(float(speedup))
        cpu_time_array.append(float(cpu_time))
        gpu_time_array.append(float(gpu_time))


    print(N, dimx, dimy, np.mean(speedup_array))
    results.append([N, dimx, dimy, np.mean(speedup_array)])


print('')
with open('{}/results_file'.format(cwd), 'w') as f:
    f.write('N, blockX, blockY, avg_time\n')
    results = sorted(results)
    for res in results:
        f.write('{}, \t{}, \t{}, \t{}\n'.format(res[0], res[1], res[2], res[3]))
