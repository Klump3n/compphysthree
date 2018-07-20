#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def computeR(MagL,MagLErr, Mag2L,Mag2LErr):
    R = Mag2L/MagL#np.divide(Mag2L,MagL)
    tmp1 = np.square(Mag2LErr/MagL)#np.square(np.divide(Mag2LErr,MagL))
    tmp2 = np.square(MagLErr*Mag2L/np.square(MagL))#np.square(np.divide(np.multiply(MagLErr,Mag2L),np.square(Mag2L)))
    RErr = np.sqrt(tmp1+tmp2)
    return R, RErr

kappa = np.arange(0.5,0.825,0.025)

#L1=16, 2L1=32, 2dim
MagL1 = np.array([0.667305,
0.794726,
0.871699,
0.980605,
1.066333,
1.145539,
1.216232,
1.291967,
1.381059,
1.446640,
1.516799,
1.587970,
1.648937
                ])

Mag2L1 = np.array([
0.571412,
0.703570,
0.804034,
0.871631,
0.990307,
1.068324,
1.133707,
1.231451,
1.286406,
1.378998,
1.440305,
1.522194,
1.573214
                ])

MagL1Err = np.array([0.015789,
0.007096,
0.011382,
0.008205,
0.009027,
0.006004,
0.009243,
0.009848,
0.006367,
0.006111,
0.005176,
0.004770,
0.005805
                     ])


Mag2L1Err = np.array([0.013786,
0.009052,
0.012185,
0.026703,
0.010321,
0.015221,
0.018902,
0.014714,
0.016531,
0.012907,
0.011807,
0.008789,
0.007334
                      ])

R1,R1Err = computeR(MagL1,MagL1Err,Mag2L1,Mag2L1Err)



#L2=20 , 2L2=40, 2dim
MagL2 = np.array([0.642329,
0.758397,
0.851562,
0.937912,
1.049806,
1.114566,
1.200576,
1.275266,
1.339198,
1.418474,
1.494576,
1.552690,
1.616621
                  ])
Mag2L2 = np.array([0.504137,
0.662260,
0.779596,
0.872630,
0.948652,
1.052062,
1.115023,
1.183632,
1.282305,
1.325360,
1.379387,
1.465470,
1.551686
                     ])
MagL2Err = np.array([0.014798,
0.009487,
0.007759,
0.012951,
0.008859,
0.010653,
0.007731,
0.006321,
0.009600,
0.007215,
0.006775,
0.006576,
0.006632
                    ])


Mag2L2Err = np.array([0.021386,
0.013560,
0.012611,
0.010382,
0.014333,
0.011706,
0.013971,
0.009359,
0.008856,
0.018712,
0.016437,
0.017025,
0.006856
                    ])


R2,R2Err = computeR(MagL2,MagL2Err,Mag2L2,Mag2L2Err)

#L3=32 , 2L3=64, 2dim
MagL3 = np.array([0.533904,
0.682706,
0.804611,
0.906373,
0.990946,
1.073077,
1.160615,
1.220157,
1.297023,
1.370425,
1.444426,
1.518923,
1.565635
                  ])
Mag2L3 = np.array([0.486905,
0.362818,
0.708922,
0.793185,
0.901318,
0.588909,
1.083170,
1.157518,
1.251010,
1.276249,
0.594393,
1.448804,
1.487048
                     ])
MagL3Err = np.array([0.013800,
0.019430,
0.007737,
0.010352,
0.014520,
0.008325,
0.005789,
0.012381,
0.013445,
0.008342,
0.012382,
0.010691,
0.008631
                    ])
Mag2L3Err = np.array([0.010648,
0.100608,
0.021079,
0.022025,
0.014505,
0.183838,
0.024850,
0.011543,
0.005209,
0.021829,
0.054269,
0.011382,
0.014690
                    ])
R3,R3Err = computeR(MagL3,MagL3Err,Mag2L3,Mag2L3Err)

R_ana = 2**(-1/(4*3.141592*kappa))


plt.figure()
plt.errorbar(kappa ,R1,R1Err,fmt='bo',barsabove=True,marker='o',label='L=16')
plt.errorbar(kappa,R2, R2Err,fmt='gx',marker='x',barsabove=True,label='L=20')
plt.errorbar(kappa,R3, R3Err,fmt='yd',marker='d',barsabove=True,label='L=32')
plt.plot(kappa, R_ana,'r',label='analytisch')
plt.title('Konvergenzverhalten von $R(L)$')
plt.xlabel('$\kappa$')
plt.ylabel('R')
plt.legend(loc='best')
plt.tight_layout()
#plt.yscale('log')
#plt.grid(1)
plt.savefig('../figures/conv_2d.pdf')
