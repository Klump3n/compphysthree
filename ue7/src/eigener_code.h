#ifndef EIGENER_CODE
#define EIGENER_CODE

cuDoubleComplex mag(cuDoubleComplex *phi);
void spin_update(double delta, double lambda, double kappa, cuDoubleComplex h);
cuDoubleComplex spin_proposal(int idx, double delta, double rand1, double rand2);
double p_a(int idx, double delta, double lambda, double kappa, cuDoubleComplex h, double rand1, double rand2);
double delta_fitting_step(int fit_count, double delta, double lambda, double kappa, cuDoubleComplex h);
double delta_fitting(double delta, double lambda, double kappa, cuDoubleComplex h);
bool spin_update_one_point(
                           int idx, double delta,
                           double lambda, double kappa, cuDoubleComplex h,
                           double rand1, double rand2, double rand3
                           );

#endif
