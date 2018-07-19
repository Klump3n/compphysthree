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
#L=4, 2L=8, 4dim
MagL1 = np.array([
                 0.0032004,
                 0.0036077,
                 0.0045948,
                 0.0057298,
                 0.0078154,
                 0.0114242,
                 0.0271091,
                 0.1066282,
                 0.4394290,
                 0.6946020,
                 0.9042737,
                 1.0831917,
                 1.2282713,
                 1.3690923,
                 1.504678
                ])

Mag2L1 = np.array([
                    0.000176,
                    0.000222,
                    0.000267,
                    0.000353,
                    0.000479,
                    0.000813,
                    0.001509,
                    0.010995,
                    0.412499,
                    0.680849,
                    0.888492,
                    1.062425,
                    1.217979,
                    1.363627,
                    1.498702
                ])

MagL1Err = np.array([
                 0.000134,
                 0.000157,
                 0.000218,
                 0.000278,
                 0.000454,
                 0.000732,
                 0.002281,
                 0.011832,
                 0.008160,
                 0.007010,
                 0.004617,
                 0.003917,
                 0.003733,
                 0.003423,
                 0.003245
                     ])


Mag2L1Err = np.array([
                     0.000008,
                     0.000010,
                     0.000011,
                     0.000019,
                     0.000026,
                     0.000058,
                     0.000133,
                     0.001848,
                     0.001885,
                     0.001412,
                     0.001090,
                     0.001027,
                     0.000916,
                     0.001022,
                     0.000914
                      ])

R1,R1Err = computeR(MagL1,MagL1Err,Mag2L1,Mag2L1Err)



#L=6 , 2L=12, 4dim
MagL2 = np.array([
                0.000615,
                0.000654,
                0.000890,
                0.001080,
                0.001578,
                0.002226,
                0.004894,
                0.035974,
                0.419866,
                0.686194,
                0.889607,
                1.066808,
                1.222554,
                1.366453,
                1.503673
                  ])
Mag2L2 = np.array([
                     0.000035,
                     0.000046,
                     0.000053,
                     0.000066,
                     0.000090,
                     0.000144,
                     0.000317,
                     0.003333,
                     0.405924,
                     0.677640,
                     0.886014,
                     1.062961,
                     1.218784,
                     1.361856,
                     1.497598
                     ])
MagL2Err = np.array([
                     0.000028,
                     0.000029,
                     0.000041,
                     0.000055,
                     0.000090,
                     0.000168,
                     0.000545,
                     0.005009,
                     0.003921,
                     0.002681,
                     0.002520,
                     0.002055,
                     0.001550,
                     0.001473,
                     0.001401
                    ])


Mag2L2Err = np.array([
                    0.000001,
                    0.000002,
                    0.000003,
                    0.000003,
                    0.000006,
                    0.000009,
                    0.000036,
                    0.000824,
                    0.000989,
                    0.000797,
                    0.000783,
                    0.000517,
                    0.000574,
                    0.000402,
                    0.000434
                    ])


R2,R2Err = computeR(MagL2,MagL2Err,Mag2L2,Mag2L2Err)
plt.figure()
plt.errorbar(kappa,R1,R1Err,marker='o',label='L=10')
plt.errorbar(kappa, R2, R2Err,marker='x',label='L=8')
plt.title('R im 4D-System')
plt.xlabel('$\kappa$')
plt.ylabel('R')
plt.legend(loc='best')
plt.tight_layout()
plt.grid(1)
plt.savefig('../figures/4dR.pdf')
