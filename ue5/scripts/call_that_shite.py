#!/usr/bin/env python3
"""
Get metrics.

"""
import os
from subprocess import call

N = 32
dimx = 256
dimy = 1

cwd = os.getcwd()
call(['{}/../bin/run-cg'.format(cwd) ,'{}'.format(N), '{}'.format(dimx), '{}'.format(dimy)])
