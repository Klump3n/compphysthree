#include <stdlib.h>
#include <stdio.h>
#include "global.h"
/*
   ueber global.h werden folgende globale Variablen eingebunden: 
   int   ndim, nvol, *lsize, **nn;
*/

void geom_pbc(){
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
  int   i, k;
  int *ibase, *ix;

  ibase = (int *) malloc((ndim+2) * sizeof(int));
  ix = (int *) malloc((ndim+1) * sizeof(int));

  /* Basis fuer Punktindizes */
  ibase[1] = 1;
  for (k=1; k<=ndim; k++) ibase[k+1] = ibase[k]*lsize[k];
  nvol = ibase[ndim+1];

  if (nn) free(nn[0]);
  free(nn);
  nn = (int **) malloc((2*ndim+1) * sizeof(int *));
  nn[0] = (int *) malloc(nvol*(2*ndim+1) * sizeof(int));
  for (k=1; k<=2*ndim; k++) nn[k] = nn[0] + nvol*k;

  for (k=1; k<=ndim; k++) ix[k] = 0;   /* Koord. des Anfangspunkts */

  for (i=0; i<nvol; i++){           /* Schleife ueber Punkte */
    for (k=1; k<=ndim; k++){
      nn[k][i] = i + ibase[k];      /* Nachbar x + e_k */
      if (ix[k] == (lsize[k]-1)) nn[k][i] -= ibase[k+1];

      nn[ndim+k][i] = i - ibase[k]; /* Nachbar x - e_k */
      if (ix[k] == 0) nn[ndim+k][i] += ibase[k+1];
    }

    for (k=1; k<=ndim; k++){        /* Koord. des naechsten Punkts */
      ix[k]++;
      if (ix[k] < lsize[k]) break;
      ix[k] = 0;
    }
  }
  free(ibase); free(ix);
}
