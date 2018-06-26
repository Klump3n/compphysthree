#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include <assert.h>

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

int other_test() {
  int i;
  int idx = 0;
  double delta = 1e-3;
  double rand1;
  double rand2;
  for (i=0;i<10;i++){
    rand1 = (double) (rand() / 30000000000.0);
    rand2 = (double) (rand() / 30000000000.0);
    double bla = akzeptanz(idx, delta, rand1, rand2);
    printf("%e, %e, %f\n", rand1, rand2, bla);
  }

  return 0;
}
