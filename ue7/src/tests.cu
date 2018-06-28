#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include <assert.h>
#include "randgpu.h"

#include "eigener_code.h"
#include "action.h"



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

int spin_update_test() {
	
  cuDoubleComplex *phi_backup;
  phi_backup = (cuDoubleComplex *) malloc(nvol*sizeof(cuDoubleComplex));
  memcpy(phi_backup, phi, nvol*sizeof(cuDoubleComplex));
  cuDoubleComplex z = make_cuDoubleComplex(0.0, 0.0);
  int idx;

  phi[0]=make_cuDoubleComplex(0.0, 0.0);
  for (idx=1; idx<nvol; idx++){
    phi[idx] = z;
  }
   
  cuDoubleComplex h=make_cuDoubleComplex(0.0,0.0);
  double lambda=0.0;
  double kappa=0.06;
  double delta = 1e-2;

  double rand1=0.8;
  double rand2=0.8;

  double phi_abs = rand1*rand1 +rand2*rand2;
  double p_analytic = exp(-phi_abs - lambda*(phi_abs-1)*(phi_abs-1)+lambda);
  double rand3 = p_analytic;
  printf("p_analytic: %f\n",rand3);
bool res1 =  spin_update_one_point(
                              0, delta,
                             lambda, kappa, h,
                             rand1, rand2, rand3+0.01
                             );

bool res2 =  spin_update_one_point(
                              0, delta,
                             lambda, kappa, h,
                             rand1, rand2, rand3-0.01
                             );

  printf("res+:%d (=1?), res-:%d (=0?)\n",res1,res2);
  memcpy(phi, phi_backup, nvol*sizeof(cuDoubleComplex));
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
/*
int delta_fitting_test() {
  cuDoubleComplex h=make_cuDoubleComplex(0.3,0.5);
  double lambda=0.7;
  double kappa=0.06;
  double delta = 1e-5;

  double delta_result = delta_fitting(delta, lambda, kappa, h);
  printf("delta result = %f\n", delta_result);

  return 0;
}
*/
int spin_set_test() {
  cuDoubleComplex h=make_cuDoubleComplex(0.3,0.5);
  double lambda=0.7;
  double kappa=0.06;
  double delta = 1e-5;

  spin_update(delta, lambda, kappa, h);

  return 0;



}

int boltzmag() {
	
  cuDoubleComplex h=make_cuDoubleComplex(0.5,0.0);
  double lambda=1.0;
  double kappa=0.1;
  double delta = 1e-5;
  random_cnfg();

  boltzmann_exp(delta, lambda, kappa, h);
  printf("\n--------------------------------\n");
  random_cnfg();
  h=make_cuDoubleComplex(0.0,0.0);
  boltzmann_exp(delta, lambda, kappa, h);

  printf("\n--------------------------------\n");

  h=make_cuDoubleComplex(0.5,0.0);

  cuDoubleComplex z=make_cuDoubleComplex(0.5, 0.0);
  for (int idx=0; idx<nvol; idx++){
    phi[idx] = z;
  }
  
  boltzmann_exp(delta, lambda, kappa, h);

  printf("\n--------------------------------\n");
  for (int idx=0; idx<nvol; idx++){
    phi[idx] = z;
  }
  h=make_cuDoubleComplex(0.0,0.0);
  boltzmann_exp(delta, lambda, kappa, h);
  return 0;


}
