#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def computeR(MagL,MagLErr, Mag2L,Mag2LErr):
    R = Mag2L/MagL#np.divide(Mag2L,MagL)
    tmp1 = np.square(Mag2LErr/MagL)#np.square(np.divide(Mag2LErr,MagL))
    tmp2 = np.square(MagLErr*Mag2L/np.square(MagL))#np.square(np.divide(np.multiply(MagLErr,Mag2L),np.square(Mag2L)))
    RErr = np.sqrt(tmp1+tmp2)
    return R, RErr


#L1=4 , 2L=8 #######################################################################
MagL1 = np.array([ 0.014128,
                   0.004398,
                   0.037572
                  ])
Mag2L1 = np.array([
                     0.002123,
                     0.000283,
                     0.001271,
                     ])
MagL1Err = np.array([
                    0.000604,
                    0.000208,
                    0.001886
                    ])

Mag2L1Err = np.array([0.000094,
                      0.000013,
                      0.000057
                    ])
R1,R1Err = computeR(MagL1,MagL1Err,Mag2L1,Mag2L1Err)
#L2=6 , 2L2=12#######################################################################
MagL2 = np.array([ 0.014128,
                   0.004398
                  ])
Mag2L2 = np.array([
                     0.002123,
                     0.000283
                     ])
MagL2Err = np.array([0.000210,0.000040])


Mag2L2Err = np.array([0.000026,0.000003
                    ])

R2,R2Err = computeR(MagL2,MagL2Err,Mag2L2,Mag2L2Err)
#L3=8 , 2L3=16#######################################################################
MagL3 = np.array([
                    0.002014
                  ])
Mag2L3 = np.array([
                    0.000227
                     ])
MagL3Err = np.array([0.000097])

Mag2L3Err = np.array([0.000010
                    ])


R3,R3Err = computeR(MagL3,MagL3Err,Mag2L3,Mag2L3Err)

d_fine = np.arange(3,5.01,0.01)
R_ana = 2**(-1*d_fine)

plt.figure()
plt.errorbar(np.array([3,4,5]),R1,R1Err,fmt='bo',barsabove=True,marker='o',label='L=4')
plt.errorbar(np.array([3,4]), R2, R2Err,fmt='gx',marker='x',barsabove=True,label='L=6')
plt.errorbar(np.array([3]), R3, R3Err,fmt='yd',marker='d',barsabove=True,label='L=8')
plt.plot(d_fine, R_ana,'r',label='analytisch')
plt.title('Konvergenzverhalten von $R(L)$')
plt.xlabel('Dimension $d$')
plt.ylabel('R')
plt.legend(loc='best')
plt.tight_layout()
plt.yscale('log')
#plt.grid(1)
plt.savefig('../figures/conv_nd.pdf')
