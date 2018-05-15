#!/usr/bin/env python3
"""
Print the execution time for the configuration layout.

"""
import numpy as np

results = np.genfromtxt('grid_parameter_results', delimiter=',', usecols=0, dtype=float)
parameters = np.genfromtxt('grid_parameters', delimiter=',', usecols=(0, 1, 2, 3), dtype=int)

# first 40 indices
low_indices = results.argsort()[:40]

for index in low_indices:
    print("{} seconds for \t gridX: {} \t gridY: {} \t threadX: {} \t threadY: {}".format(results[index], parameters[index][0], parameters[index][1], parameters[index][2], parameters[index][3]))
