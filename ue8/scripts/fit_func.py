#!/usr/bin/env python3
"""
Fit a function.

"""
import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt


def func(kappa, factor, kappa_c, c):
    return factor * np.power((kappa - kappa_c), c)

x = [
    0.0,
    0.023329999999999997,
    0.04665999999999999,
    0.06998999999999998,
    0.0933199999999999,
    0.11665,
    0.13998,
    0.16330999999999998,
    0.18664000000000003,
    0.20997000000000002,
    0.2333,
    0.25663,
    0.27996,
    0.30329,
    0.32661999999999997,
    0.34995,
    0.37328000000000006,
    0.3966100000000001,
    0.41994000000000004,
    0.44327,
    0.4666
]

res = [
    0.012174, 0.000284,
    0.014168, 0.000307,
    0.015852, 0.00034,
    0.018891, 0.000413,
    0.021366, 0.000507,
    0.026784, 0.000687,
    0.033721, 0.000885,
    0.042882, 0.001303,
    0.057489, 0.002974,
    0.092706, 0.003319,
    0.174315, 0.006588,
    0.325252, 0.008191,
    0.531875, 0.008633,
    0.709577, 0.007085,
    0.855642, 0.006863,
    0.993551, 0.006142,
    1.098343, 0.006575,
    1.206233, 0.009107,
    1.303047, 0.011901,
    1.404452, 0.009953,
    1.497503, 0.0089
]

datalen = int(len(res)/2)

kappa = []
mean = []
std = []

for k in range(datalen):
    m = res[2*k]
    s = res[2*k + 1]

    if (m > .1):
        kappa.append(x[k])
        mean.append(m)
        std.append(s)

popt, pcov = so.curve_fit(func, kappa, mean, p0=[3, .2, .5], sigma=std)
print(popt)
print(pcov)
print(np.sqrt(pcov[0,0]))
print(np.sqrt(pcov[1,1]))
print(np.sqrt(pcov[2,2]))
plt.figure()
ebar,_,_ = plt.errorbar(kappa, mean, yerr=std, marker='_', ls='None')
print(ebar.get_color())
plt.plot(kappa, func(kappa, *popt), label=r'L = {}, '.format(4) + r'$\kappa_{c}$ = ' + '{:.4f}'.format(popt[1]), color=ebar.get_color())
plt.legend(loc='best')
plt.show()
