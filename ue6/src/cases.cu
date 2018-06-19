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
  ndim = 2;                    /* dim des problems */
  lsize = (int*) malloc((ndim + 1) * sizeof(int));
  for (int k=0; k<ndim; k++)
    {
      lsize[k+1]=10;
    };

  geom_pbc();                  /* berechne indizes nn[k][i] und volumen */

  phi=(cuDoubleComplex*)malloc(nvol*sizeof(cuDoubleComplex));
}

void aufg_1() {
  printf("Aufgabe 1\n");

  cuDoubleComplex z = make_cuDoubleComplex(2.125,1.0);
  for (int k=0; k<nvol; k++)
    {
      phi[k]=z;
    }

  cuDoubleComplex h = make_cuDoubleComplex(1.0, 1.0);
  double kappa = 1.0;
  double lambda = 1.0;

  printf("dimension des problems %d\n", nvol);

  double S_MC = S(phi, h, kappa, lambda);
  double S_ANA = S_analytical(z, h, kappa, lambda);
  assert(abs((S_MC-S_ANA)/S_ANA)<sqrt(nvol)*DBL_EPSILON);
  printf("wirkung monte carlo %f\n", S_MC);
  printf("wirkung analytisch %f\n", S_ANA);

  /* printf("p(phi(x)) %f\n", p(0, h, kappa, lambda)); */

}

void aufg_2() {
  printf("Aufgabe 2\n");

  for (int k=0; k<nvol; k++)
    {
      phi[k]= make_cuDoubleComplex(
                                   (float)(rand() & 0xFF ) / 10.0,
                                   (float)(rand() & 0xFF ) / 10.0
                                   );
    };

  cuDoubleComplex phi_x;
  double alpha = (float)(rand() & 0xFF ) / 10.0;
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
  double kappa = 1.0;
  double lambda = 1.0;

  double S_mc_orig = S(phi, h, kappa, lambda);
  double S_mc_mod = S(phi_mod, h, kappa, lambda);
  assert(abs((S_mc_orig-S_mc_mod)/S_mc_orig)<sqrt(nvol)*DBL_EPSILON);
  printf("wirkung original monte carlo %f\n", S_mc_orig);
  printf("wirkung modifiziert monte carlo %f\n", S_mc_mod);

}

void aufg_3() {
  cuDoubleComplex *phi_one = (cuDoubleComplex*)malloc(nvol*sizeof(cuDoubleComplex));
  cuDoubleComplex *phi_two = (cuDoubleComplex*)malloc(nvol*sizeof(cuDoubleComplex));
  for (int k=0; k<nvol; k++)
    {
      phi_one[k] = make_cuDoubleComplex(
                                   (float)(rand() & 0xFF ) / 100.0,
                                   (float)(rand() & 0xFF ) / 100.0
                                   );
      phi_two[k] = phi_one[k];
    };

  phi_two[1] = make_cuDoubleComplex( /* Because. */
                                                (float)(rand() & 0xFF ) / 100.0,
                                                (float)(rand() & 0xFF ) / 100.0
                                                );

  cuDoubleComplex h = make_cuDoubleComplex(1.0, 1.0);
  double kappa = 1.0;
  double lambda = 1.0;

  double S_one = S(phi_one, h, kappa, lambda);
  double S_two = S(phi_two, h, kappa, lambda);
  double P_div_S1_by_S2 = exp(S_two - S_one);

  double P_one = p(1, phi_one[1], h, kappa, lambda);
  double P_two = p(1, phi_two[1], h, kappa, lambda);

  printf("p(1) %f\n", P_one);
  printf("p(2) %f\n", P_two);
  printf("p(1)/p(2) %f\n", P_one/P_two);
  printf("P(1)/P(2) %f\n", P_div_S1_by_S2);


}
