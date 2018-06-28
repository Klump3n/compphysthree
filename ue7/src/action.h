#ifndef ACTION_H
#define ACTION_H

double alocal(int idx, double lambda, double kappa, cuDoubleComplex h);
double action(double lambda, double kappa, cuDoubleComplex h);
void random_cnfg(void);

#endif
