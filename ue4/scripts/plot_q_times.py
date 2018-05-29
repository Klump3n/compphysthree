#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('../measurements/hand_pasted_times_aufg_5', skip_header=2, delimiter=',', usecols=(0,1,2))

plt.figure()
plt.plot(data[:,0], data[:, 1], label=r'int')
plt.plot(data[:,0], data[:, 2], label=r'double')
plt.title('Vergleich der Laufzeiten bei unterschiedlichen Datentypen\nqUnrollReducing')
plt.xlabel('q-Wert')
plt.ylabel('Laufzeit der Funktion [ms]')
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('../figures/laufzeiten_qunroll_int_double.pdf')
