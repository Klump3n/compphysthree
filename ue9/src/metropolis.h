#ifndef METROPOLIS_H
#define METROPOLIS_H

extern double metro_sweep(double delta);
extern double metro_sweep_alt(double delta, int *evenArray, int *oddArray, double *rnd);
extern int update_spin(int idx, double* rnd, double del, int ntrial);

#endif
