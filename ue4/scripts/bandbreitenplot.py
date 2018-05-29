#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

unroll = np.genfromtxt('../measurements/hand_paste_2', skip_header=2, skip_footer=7, delimiter=',', usecols=(0,1,2,3,4))
interleaved = np.genfromtxt('../measurements/hand_paste_2', skip_header=11, delimiter=',', usecols=(0,1,2,3,4))

print(interleaved)
print(unroll)


index = 1
plt.figure()
plt.plot(unroll[:,0], unroll[:, index], label=r'Unroll')
plt.plot(unroll[:,0], interleaved[:, index], label=r'Interleaved')
plt.title('Vergleich der Laufzeiten zwischen\n unroll und interleaved')
plt.xlabel('Blocksize')
plt.ylabel('Laufzeit der Funktion [ms]')
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('../figures/laufzeiten_blocksize.pdf')

index = 2
plt.figure()
plt.plot(unroll[:,0], unroll[:, index], label=r'Unroll')
plt.plot(unroll[:,0], interleaved[:, index], label=r'Interleaved')
plt.title('Vergleich der memory load efficiency zwischen\n unroll und interleaved')
plt.xlabel('Blocksize')
plt.ylabel('Memory Load Efficiency [%]')
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('../figures/mem_ld_eff_blocksize.pdf')

index = 3
plt.figure()
plt.plot(unroll[:,0], unroll[:, index], label=r'Unroll')
plt.plot(unroll[:,0], interleaved[:, index], label=r'Interleaved')
plt.title('Vergleich der Bandbreite zwischen\n unroll und interleaved')
plt.xlabel('Blocksize')
plt.ylabel('Bandbreite [GB/s]')
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('../figures/bw_blocksize.pdf')

index = 4
plt.figure()
plt.plot(unroll[:,0], unroll[:, index]*100, label=r'Unroll')
plt.plot(unroll[:,0], interleaved[:, index]*100, label=r'Interleaved')
plt.title('Vergleich der Auslastung zwischen\n unroll und interleaved')
plt.xlabel('Blocksize')
plt.ylabel('Auslastung [%]')
plt.xscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('../figures/occ_blocksize.pdf')
