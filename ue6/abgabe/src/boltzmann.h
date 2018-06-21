#ifndef BOLTZMANN_H
#define BOLTZMANN_H

double S(cuDoubleComplex *phi, cuDoubleComplex h, double kappa, double lambda);
double S_analytical(cuDoubleComplex z, cuDoubleComplex h, double kappa, double lambda);
double p(int x, cuDoubleComplex *phi, cuDoubleComplex h, double kappa, double lambda);
double p_arg(int x, cuDoubleComplex *phi, cuDoubleComplex h, double kappa, double lambda);

#endif
