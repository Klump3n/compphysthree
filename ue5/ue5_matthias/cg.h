#ifndef CG_H
#define CG_H

extern double cg(double *x, double *r, int maxiter, double rel, int *status);

extern double cg_gpu(double *x, double *r, int maxiter, double rel, int *status);

#endif
