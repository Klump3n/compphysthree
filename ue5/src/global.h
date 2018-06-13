

#ifndef GLOBAL_H
#define GLOBAL_H

#define DBL_EPSILON 2.2204460492503131e-16

#ifdef MAIN_PROGRAM
   #define EXTERN
#else
   #define EXTERN extern
#endif

/*
   Globale Variablen stehen in allen Funktionen zur Verfuegung.
   Achtung: Das gilt *nicht* fuer Kernel-Funktionen!
*/
EXTERN int Nx, Ny, npts;
EXTERN int *active;
EXTERN dim3 block, grid;
EXTERN int nblk;
EXTERN dim3 block2, grid2;      /* unroll stage 1 */
EXTERN dim3 block3, grid3;      /* unroll stage 2 */

#undef EXTERN

#endif
