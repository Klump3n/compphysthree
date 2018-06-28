#ifndef EIGENER_CODE
#define EIGENER_CODE

cuDoubleComplex mag(cuDoubleComplex *phi);
void spin_update(double delta, double lambda, double kappa, cuDoubleComplex h);
cuDoubleComplex spin_proposal(int idx, double delta, double rand1, double rand2);
double p_a(int idx, double delta, double lambda, double kappa, cuDoubleComplex h, double rand1, double rand2);
bool spin_update_one_point(
                           int idx, double delta,
                           double lambda, double kappa, cuDoubleComplex h,
                           double rand1, double rand2, double rand3
                           );
void boltzmann_exp(double delta, double lambda, double kappa, cuDoubleComplex h); 

#endif
