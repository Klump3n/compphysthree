#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def computeR(MagL,MagLErr, Mag2L,Mag2LErr):
    R = Mag2L/MagL#np.divide(Mag2L,MagL)
    tmp1 = np.square(Mag2LErr/MagL)#np.square(np.divide(Mag2LErr,MagL))
    tmp2 = np.square(MagLErr*Mag2L/np.square(MagL))#np.square(np.divide(np.multiply(MagLErr,Mag2L),np.square(Mag2L)))
    RErr = np.sqrt(tmp1+tmp2)
    return R, RErr

kappa = np.array([0.00, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2,
                    0.225, 0.25, 0.275, 0.3, 0.325, 0.35])
#L=10, 2L=20, 3dim

MagL1 = np.array([
                 0.000781,
                 0.000959,
                 0.001122,
                 0.001182,
                 0.001406,
                 0.001733,
                 0.002361,
                 0.003289,
                 0.006624,
                 0.012158,
                 0.067340,
                 0.405582,
                 0.614301,
                 0.783955,
                 0.941954
                ])

Mag2L1 = np.array([
                     0.000092,
                     0.000110,
                     0.000122,
                     0.000145,
                     0.000185,
                     0.000211,
                     0.000314,
                     0.000422,
                     0.000541,
                     0.001608,
                     0.015514,
                     0.372903,
                     0.600169,
                     0.776499,
                     0.919605
                ])

MagL1Err = np.array([
                     0.000034,
                     0.000049,
                     0.000060,
                     0.000057,
                     0.000072,
                     0.000091,
                     0.000143,
                     0.000218,
                     0.000556,
                     0.001718,
                     0.010434,
                     0.006675,
                     0.004660,
                     0.004760,
                     0.003026
                     ])


Mag2L1Err = np.array([
                    0.000004,
                    0.000005,
                    0.000005,
                    0.000007,
                    0.000009,
                    0.000012,
                    0.000018,
                    0.000032,
                    0.000052,
                    0.000356,
                    0.004838,
                    0.003478,
                    0.001873,
                    0.002710,
                    0.002473
                      ])

R1,R1Err = computeR(MagL1,MagL1Err,Mag2L1,Mag2L1Err)



#L=8 , 2L=16, 3dim
MagL2 = np.array([
                 0.001398,
                 0.001802,
                 0.002044,
                 0.002278,
                 0.002527,
                 0.003013,
                 0.004350,
                 0.006443,
                 0.009740,
                 0.019025,
                 0.117825,
                 0.395040,
                 0.621997,
                 0.794447,
                 0.942183
                  ])
Mag2L2 = np.array([
                    0.000188,
                    0.000219,
                    0.000251,
                    0.000283,
                    0.000385,
                    0.000412,
                    0.000567,
                    0.000699,
                    0.001288,
                    0.002686,
                    0.037330,
                    0.373866,
                    0.597733,
                    0.770346,
                    0.922988
                                     ])
MagL2Err = np.array([
                     0.000065 ,
                     0.000086,
                     0.000090,
                     0.000110,
                     0.000144,
                     0.000155,
                     0.000313,
                     0.000472,
                     0.000693,
                     0.002148,
                     0.011558,
                     0.009935,
                     0.005298,
                     0.003842,
                     0.004178
                    ])


Mag2L2Err = np.array([
                     0.000008,
                     0.000010,
                     0.000012,
                     0.000014,
                     0.000019,
                     0.000024,
                     0.000030,
                     0.000050,
                     0.000124,
                     0.000304,
                     0.009071,
                     0.004040,
                     0.002275,
                     0.003539,
                     0.001913
                      ])


R2,R2Err = computeR(MagL2,MagL2Err,Mag2L2,Mag2L2Err)
plt.figure()
plt.errorbar(kappa,R1,R1Err,marker='o',label='L=10')
plt.errorbar(kappa, R2, R2Err,marker='x',label='L=8')
plt.title('R im 3D-System')
plt.xlabel('$\kappa$')
plt.ylabel('R')
plt.legend(loc='best')
plt.tight_layout()
plt.grid(1)
plt.savefig('../figures/3dR.pdf')
