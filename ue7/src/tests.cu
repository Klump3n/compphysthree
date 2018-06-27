#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include <assert.h>
#include "randgpu.h"

#include "eigener_code.h"


int mag_test() {

  cuDoubleComplex *phi_test;
  phi_test = (cuDoubleComplex *) malloc(nvol*sizeof(cuDoubleComplex));
  cuDoubleComplex z = make_cuDoubleComplex(1.0, 2.0);
  int idx;

  for (idx=0; idx<nvol; idx++){
    phi_test[idx] = z;
  }

  cuDoubleComplex phi_result = mag(phi_test);

  assert(cuCabs(cuCsub(phi_result, z))/cuCabs(z)<nvol*DBL_EPSILON);
  printf("Check magnetisierung erfolgreich.\n");
  printf("%e = |M - z| / |z| < nvol*DBL_EPSILON = %e\n", cuCabs(cuCsub(phi_result, z))/cuCabs(z), nvol*DBL_EPSILON);
  return 0;
}

int update_test() {
  /* generate random numbers */
  return 0;
}

int other_test() {

  cuDoubleComplex h=make_cuDoubleComplex(0.3,0.5);
  double lambda=0.7;
  double kappa=0.06;
  double delta = 1e-1;

  double *random_nums;
  random_nums = (double *) malloc(2*sizeof(double));
  random_nums = randgpu(2);
  double rand1 = random_nums[0],
    rand2 = random_nums[1];

  int idx = 0;

  double p_a_val = p_a(idx, delta, lambda, kappa, h, rand1, rand2);
  printf("%e, %e, %f\n", rand1, rand2, p_a_val);

  return 0;
}

int delta_fitting_test() {
  cuDoubleComplex h=make_cuDoubleComplex(0.3,0.5);
  double lambda=0.7;
  double kappa=0.06;
  double delta = 1e-5;

  double delta_result = delta_fitting(delta, lambda, kappa, h);
  printf("delta result = %f\n", delta_result);

  return 0;
}

int spin_set_test() {
  cuDoubleComplex h=make_cuDoubleComplex(0.3,0.5);
  double lambda=0.7;
  double kappa=0.06;
  double delta = 1e-5;

  spin_update(delta, lambda, kappa, h);

  return 0;

}
