#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include "global.h"

double S(cuDoubleComplex *phi, cuDoubleComplex h, double kappa, double lambda) {

  /*
    Folgende Groessen muessen gesetzt sein:             B Bunk 12/2005
    Dimension     ndim                              rev     4/2013
    Gittergroesse lsize[k], k=1..ndim

    Angelegt und berechnet wird
    Volumen       nvol
    NN-Indexfeld  nn[k][i], k=0..2*ndim, i=0..(nvol-1)

    nn[k][i] gibt den Index des Nachbarn von i in Richtung +k,
    unter Beruecksichtigung periodischer Randbedingungen.
    Fuer einen Schritt in Richtung -k setze man den Index auf (ndim+k).
    nn[0][i] ist reserviert.
  */
  double S_val = 0.0;

  double kappa_sum = 0.0;
  double phi_norm;

  /* nvol aus global */
  for (int i=0; i<nvol; i++)
    {

      kappa_sum = 0.0;

      for (int j=1; j<(ndim+1); j++) {
        kappa_sum += 2 * (cuCreal(phi[i]) * cuCreal(phi[ nn[j][i] ]) +
                          cuCimag(phi[i]) * cuCimag(phi[ nn[j][i] ]));
      }

      phi_norm = cuCabs(phi[i])*cuCabs(phi[i]);

      S_val +=
        phi_norm +
        lambda * (phi_norm - 1) * (phi_norm - 1) -
        kappa * kappa_sum -
        2 * (cuCreal(h) * cuCreal(phi[i]) +
             cuCimag(h) * cuCimag(phi[i]));
        }

  return S_val;
}

double S_analytical(cuDoubleComplex z, cuDoubleComplex h, double kappa, double lambda) {
  return nvol * (
                 (1 - 2*kappa*ndim) * cuCabs(z)*cuCabs(z) +
                 lambda * (cuCabs(z)*cuCabs(z) - 1) * (cuCabs(z)*cuCabs(z) - 1) -
                 2 * (cuCreal(h) * cuCreal(z) +
                      cuCimag(h) * cuCimag(z))
                 );
}

double p_arg(int x, cuDoubleComplex *phi, cuDoubleComplex h, double kappa, double lambda) {

  cuDoubleComplex phi_x = phi[x];

  double phi_norm = cuCabs(phi_x)*cuCabs(phi_x);
  cuDoubleComplex kappa_sum = make_cuDoubleComplex(0.0, 0.0);

  for (int j=1; j<(ndim + 1); j++) {
    kappa_sum = cuCadd(kappa_sum,
                       cuCadd(phi[ nn[j][x] ], phi[ nn[ndim+j][x] ])
                       );
  };
  kappa_sum = cuCmul(make_cuDoubleComplex(kappa, 0.0), kappa_sum);
  cuDoubleComplex Bx = cuCadd(h, kappa_sum);

  return 2 * (cuCreal(Bx) * cuCreal(phi_x) + cuCimag(Bx) * cuCimag(phi_x)) -
    phi_norm - lambda * (phi_norm - 1) * (phi_norm - 1);
}

double p(int x, cuDoubleComplex *phi, cuDoubleComplex h, double kappa, double lambda) {
  return exp(p_arg(x, phi, h, kappa, lambda));
}
