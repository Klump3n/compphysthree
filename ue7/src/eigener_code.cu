#include "global.h"

#include "eigener_code.h"
#include "randgpu.h"
#include "action.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>




cuDoubleComplex mag(cuDoubleComplex *phi) {
  int idx;
  cuDoubleComplex magRes = make_cuDoubleComplex(0.0, 0.0);

  for (idx=0; idx<nvol; idx++){
    magRes = cuCadd(magRes, phi[idx]);
  }

  return cuCdiv(magRes, make_cuDoubleComplex((double)nvol, 0.0));
}


void spin_update(double delta, double lambda, double kappa, cuDoubleComplex h) {
  int idx, jdx;

  int i = 0,
    max_tries = 1000;

  double *random_nums;
  /* 2 for the random vector, 1 for comparison to accept */
  /* also 10 for every sweep and nvol for every datapoint */
  random_nums = (double *) malloc(3*10*nvol*sizeof(double));

  double rand1, rand2, rand3;
  bool spin_updated = false;

  int count_success = 0;
  double accept_percentage = 0.0;

  double *rnd;

  while (((accept_percentage < 0.35) || (accept_percentage > 0.45)) && i<max_tries)
    {

      i++;
      printf("delta = %f \t", delta);
      rnd = randgpu(3*10*nvol);
      memcpy(random_nums, rnd, 3*10*nvol*sizeof(double));

      count_success = 0;
      accept_percentage = 0.0;

      /* sweep */
      for (idx=0; idx<nvol; idx++) {

        /* every point ten times */
        for (jdx=0; jdx<10; jdx++) {
          printf("AM I CORRECT?\n");
          rand1 = random_nums[3*10*idx + 3*jdx + 0]; /* I GUESS... */
          rand2 = random_nums[3*10*idx + 3*jdx + 1];
          rand3 = random_nums[3*10*idx + 3*jdx + 2];

          spin_updated = spin_update_one_point(
                                               idx, delta,
                                               lambda, kappa, h,
                                               rand1, rand2, rand3
                                               );

          if (spin_updated) {
            count_success++;
          }

          /* reset */
          spin_updated = false;

        }

      }

      printf("Counts: %d/%d\t", count_success, 10*nvol);

      accept_percentage = ((double) count_success / (double) (10*nvol));

      printf("Accepted = %f\n", accept_percentage);

      if (accept_percentage < 0.35) {
        delta = 1.05*delta;
      }
      else if (accept_percentage > 0.45) {
        delta = 0.95*delta;
      }
      else {
        break;
      }
    }

  free(random_nums);

}


double delta_fitting(double delta, double lambda, double kappa, cuDoubleComplex h)
{

  printf("DAS IST NICHT NACH AUFGABENSTELLUNG\n");

  /* how many runs to do */
  int fit_count = 10;
  fit_count = (fit_count >= nvol) ? nvol : fit_count;
  /* printf("doing %d runs\n", fit_count); */

  /* printf("\nStarting delta fitting\n"); */

  int idx;

  double *random_nums;
  /* 2 for the random vector, 1 for comparison to accept */
  random_nums = (double *) malloc(3*fit_count*sizeof(double));

  double rand1, rand2, rand3;
  bool spin_updated = false;

  int count_success = 0;
  double accept_percentage = 0.0;

  double *rnd;

  while ((accept_percentage < 0.35) || (accept_percentage > 0.45))
    {
      /* printf("Delta fitting with delta = %f\t", delta); */

      rnd = randgpu(3*fit_count);
      memcpy(random_nums, rnd, 3*fit_count*sizeof(double));

      count_success = 0;
      accept_percentage = 0.0;

      for (idx=0; idx<fit_count; idx++) {
        rand1 = random_nums[3*idx + 0];
        rand2 = random_nums[3*idx + 1];
        rand3 = random_nums[3*idx + 2];

        spin_updated = spin_update_one_point(
                                             idx, delta,
                                             lambda, kappa, h,
                                             rand1, rand2, rand3
                                             );

        if (spin_updated) {
          count_success++;
        }

        /* reset */
        spin_updated = false;

      }

      accept_percentage = ((double) count_success / (double) fit_count);

      /* printf("Accepted = %f\n", accept_percentage); */

      if (accept_percentage < 0.35) {
        delta = 1.05*delta;
      }
      else if (accept_percentage > 0.45) {
        delta = 0.95*delta;
      }
      else {
        break;
      }
    }

  free(random_nums);
  return delta;

}

/*
 * return a new local phi
 */
cuDoubleComplex spin_proposal(int idx, double delta, double rand1, double rand2) {
  double r1 = 2 * delta * (rand1 - 0.5),
    r2 = 2 * delta * (rand2 - 0.5);
  return cuCadd(phi[idx], make_cuDoubleComplex(r1, r2));
}

/*
 * calculate p_a at idx for two random numbers
 */
double p_a(int idx, double delta, double lambda, double kappa, cuDoubleComplex h, double rand1, double rand2) {

  /* calculate local action at idx */
  double aloc_old_phi = alocal(idx, lambda, kappa, h);

  /* save old phi at idx */
  cuDoubleComplex old_phi = phi[idx];

  /* set phi at idx to a proposed new state */
  phi[idx] = spin_proposal(idx, delta, rand1, rand2);

  /* calculate action with new phi */
  double aloc_new_phi = alocal(idx, lambda, kappa, h);

  /* reset old state */
  phi[idx] = old_phi;

  /* return p_a */
  if (aloc_new_phi <= aloc_old_phi) {
    return 1.0;
  } else {
    return exp(aloc_old_phi - aloc_new_phi);
  }
}

/*
 * update a local coordinate
 */
bool spin_update_one_point(
                           int idx, double delta,
                           double lambda, double kappa, cuDoubleComplex h,
                           double rand1, double rand2, double rand3
                           )
{

  printf("TEST SCHREIBEN FUER SPIN UPDATE\n");

  /* calculate p_a for the index and the proposed random numbers */
  double p_a_val = p_a(idx, delta, lambda, kappa, h, rand1, rand2);

  /*
   * compare p_a to rand3
   * if p_a <= rand3 accept
   * else reject (or do nothing)
   * since the random numbers have equal likelihood this means that we accept a
   * number with p_a_val percentage
   */
  if (p_a_val <= rand3) {
    phi[idx] = spin_proposal(idx, delta, rand1, rand2);
    return true;
  } else {
    return false;
  }
}
