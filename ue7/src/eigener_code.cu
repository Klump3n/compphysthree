#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"

#include "randgpu.h"
#include "action.h"


cuDoubleComplex mag(cuDoubleComplex *phi) {
  int idx;
  cuDoubleComplex magRes = make_cuDoubleComplex(0.0, 0.0);

  for (idx=0; idx<nvol; idx++){
    magRes = cuCadd(magRes, phi[idx]);
  }

  return cuCdiv(magRes, make_cuDoubleComplex((double)nvol, 0.0));
}

cuDoubleComplex spin_update(cuDoubleComplex *phi, double delta) {
  int idx;

  double *random_nums;
  random_nums = (double *) malloc(2*nvol*sizeof(double));
  random_nums = randgpu(2*nvol);

  for (idx=0; idx<nvol; idx++) {
    /* phi[idx] = phi[idx] + make_cuDoubleComplex( */
    /*                                            delta * (2 * random_nums[2*idx] + 1), */
    /*                                            delta * (2 * random_nums[2*idx + 1] + 1) */
    /*                                            ); */
  }

  return *phi;
}

/* cuDoubleComplex sweep() { */
  
/* } */

cuDoubleComplex spin_update_one_point(int idx, cuDoubleComplex *phi, double delta, double rand1, double rand2) {
  return cuCadd(phi[idx],
                make_cuDoubleComplex(
                                     delta * (2 * rand1 + 1),
                                     delta * (2 * rand2 + 1)
                                     )
                );
}

double akzeptanz(int idx, double delta, double rand1, double rand2) {

  cuDoubleComplex old_phi = phi[idx];

  double lambda = 1.0;
  double kappa = 1.0;
  cuDoubleComplex h = make_cuDoubleComplex(1.0, 1.0);

  double aloc_old_phi = alocal(idx, lambda, kappa, h);
  printf("aloc_old_phi = %f\n", aloc_old_phi);

  phi[idx] = spin_update_one_point(idx, phi, delta, rand1, rand2);

  double aloc_new_phi = alocal(idx, lambda, kappa, h);
  printf("aloc_new_phi = %f\n", aloc_new_phi);

  if (aloc_old_phi <= aloc_new_phi) {
    return 1.0;
  } else {
    phi[idx] = old_phi;
    return exp(aloc_new_phi - aloc_old_phi);
  }
}
