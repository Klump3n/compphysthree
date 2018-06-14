#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

Ndata = np.array([32,64,128,256,512,1024])
avg_speedup = np.array([0.0260783,0.0895114,0.3413917,0.851386,1.4301993,2.3417831])

plt.figure()
plt.plot(Ndata, avg_speedup)#, label=r'Speedup')
#plt.title('Durchschnittlicher Speedup')
plt.xlabel('Gittergroesse N')
plt.ylabel('Durchschnittlicher Speedup')
#plt.xscale('log')
#plt.legend(loc='best')
plt.tight_layout()
plt.savefig('speedup.pdf')
