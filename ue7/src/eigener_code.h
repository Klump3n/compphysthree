#ifndef EIGENER_CODE
#define EIGENER_CODE

/* #include "global.h" */

cuDoubleComplex mag(cuDoubleComplex *phi);
cuDoubleComplex spin_update(cuDoubleComplex *phi, double delta);
cuDoubleComplex spin_update_one_point(int idx, cuDoubleComplex *phi, double delta, double rand1, double rand2);
double akzeptanz(int idx, double delta, double rand1, double rand2);

#endif
