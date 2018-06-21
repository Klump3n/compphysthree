#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include "global.h"
#include <assert.h>

#include "geom_pbc.h"
#include "boltzmann.h"

void prepare_run(){
  printf("Vorbereitung\n");

  /* ist aus global importiert */
  ndim = 3;                    /* dim des problems */
  lsize = (int*) malloc((ndim + 1) * sizeof(int));
  for (int k=0; k<ndim; k++)
    {
      lsize[k+1]=5;
    };

  geom_pbc();                  /* berechne indizes nn[k][i] und volumen */

  phi=(cuDoubleComplex*)malloc(nvol*sizeof(cuDoubleComplex));

}

void aufg_1() {
  printf("\nAufgabe 1\n");

  printf("Iteration ueber 5 verschiedene Kombinationen\n");

  printf("volumen des problems %d\n", nvol);
  printf("dimension des problems %d\n", ndim);

  double rand_scaling = 100000.0;

  for (int i = 0; i<5; i++) {

    cuDoubleComplex z = make_cuDoubleComplex(
                                             (float)(rand() & 0xFF ) / rand_scaling,
                                             (float)(rand() & 0xFF ) / rand_scaling
                                             );
    for (int k=0; k<nvol; k++)
      {
        phi[k]=z;
      }

    cuDoubleComplex h = make_cuDoubleComplex(
                                             (float)(rand() & 0xFF ) / rand_scaling,
                                             (float)(rand() & 0xFF ) / rand_scaling
                                             );
    double kappa = (float)(rand() & 0xFF ) / rand_scaling;
    double lambda = (float)(rand() & 0xFF ) / rand_scaling;

    double S_MC = S(phi, h, kappa, lambda);
    double S_ANA = S_analytical(z, h, kappa, lambda);
    printf("wirkung monte carlo %f\n", S_MC);
    printf("wirkung analytisch %f\n", S_ANA);

    printf("\nassert:\n");
    printf("abs((S_MC - S_ANA) / S_ANA) = %e\n", abs((S_MC-S_ANA)/S_ANA));
    assert(abs((S_MC-S_ANA)/S_ANA)<nvol*DBL_EPSILON);

  };

}

void aufg_2() {
  printf("\nAufgabe 2\n");

  double rand_scaling = 100000.0;

  for (int k=0; k<nvol; k++)
    {
      phi[k]= make_cuDoubleComplex(
                                   (float)(rand() & 0xFF ) / rand_scaling,
                                   (float)(rand() & 0xFF ) / rand_scaling
                                   );
    };

  cuDoubleComplex phi_x;
  double alpha = (float)(rand() & 0xFF ) / rand_scaling;
  cuDoubleComplex *phi_mod = (cuDoubleComplex*)malloc(nvol*sizeof(cuDoubleComplex));
  for (int k=0; k<nvol; k++)
    {
      phi_x = phi[k];
      phi_mod[k]=cuCadd(
                    cuCmul(phi_x, make_cuDoubleComplex(cos(alpha), 0.0)),
                    cuCmul(phi_x, make_cuDoubleComplex(0.0, sin(alpha)))
                    );
    };

  cuDoubleComplex h = make_cuDoubleComplex(0.0, 0.0);
  double kappa = (float)(rand() & 0xFF ) / rand_scaling;
  double lambda = (float)(rand() & 0xFF ) / rand_scaling;

  double S_mc_orig = S(phi, h, kappa, lambda);
  double S_mc_mod = S(phi_mod, h, kappa, lambda);
  printf("originale wirkung %f\n", S_mc_orig);
  printf("modifizierte wirkung %f\n", S_mc_mod);

  printf("\nassert:\n");
  printf("abs((S_mc_orig - S_mc_mod) / S_mc_orig) = %e\n", abs((S_mc_orig-S_mc_mod)/S_mc_orig));
  assert(abs((S_mc_orig-S_mc_mod)/S_mc_orig)<nvol*DBL_EPSILON);

}

void aufg_3() {
  printf("\nAufgabe 3\n");

  int diff_index = (int)(rand() % nvol); /* random number from 0 to (nvol - 1) */
  double rand_scaling = 100000.0;

  cuDoubleComplex *phi_one = (cuDoubleComplex*)malloc(nvol*sizeof(cuDoubleComplex));
  cuDoubleComplex *phi_two = (cuDoubleComplex*)malloc(nvol*sizeof(cuDoubleComplex));
  for (int k=0; k<nvol; k++)
    {
      phi_one[k] = make_cuDoubleComplex(
                                   (float)(rand() & 0xFF ) / rand_scaling,
                                   (float)(rand() & 0xFF ) / rand_scaling
                                   );
      phi_two[k] = phi_one[k];
    };

  phi_two[diff_index] = make_cuDoubleComplex( /* Because. */
                                             (float)(rand() & 0xFF ) / rand_scaling,
                                             (float)(rand() & 0xFF ) / rand_scaling
                                              );

  cuDoubleComplex h = make_cuDoubleComplex(1.0, 1.0);
  double kappa = 1.0;
  double lambda = 1.0;

  double S_one = S(phi_one, h, kappa, lambda);
  double S_two = S(phi_two, h, kappa, lambda);
  double P_one = exp(-S_one);
  double P_two = exp(-S_two);

  double P_div_S1_by_S2 = exp(S_two - S_one);

  double p_one_arg = p_arg(diff_index, phi_one, h, kappa, lambda);
  double p_two_arg = p_arg(diff_index, phi_two, h, kappa, lambda);
  double p_one = p(diff_index, phi_one, h, kappa, lambda);
  double p_two = p(diff_index, phi_two, h, kappa, lambda);

  double div_p1_by_p2 = exp(p_one_arg - p_two_arg);

  printf("p(1) %f\n", p_one);
  printf("p(2) %f\n", p_two);

  printf("P(1) %f\n", P_one);
  printf("P(2) %f\n", P_two);

  printf("S(1) %f\n", S_one);
  printf("S(2) %f\n", S_two);

  printf("p(1)/p(2) %f\n", p_one/p_two);
  printf("P(1)/P(2) %f\n", P_one/P_two);

  printf("p(1)/p(2) alternativ %f\n", div_p1_by_p2);
  printf("P(1)/P(2) alternativ %f\n", P_div_S1_by_S2);

  printf("\nassert:\n");
  printf("abs((div_p1_by_p2 - P_div_S1_by_S2) / div_p1_by_p2) = %e\n", abs((div_p1_by_p2 - P_div_S1_by_S2)/div_p1_by_p2));
  assert(abs((div_p1_by_p2 - P_div_S1_by_S2)/div_p1_by_p2)<nvol*DBL_EPSILON);
}
