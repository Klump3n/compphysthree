#ifndef ADDED_STUFF_GPU_H
#define ADDED_STUFF_GPU_H

#include <cuComplex.h>

typedef cuDoubleComplex spin;

#define make_spin(a,b) (                          \
                        make_cuDoubleComplex(a,b) \
                         )

int gpu_sweep(spin *d_phi,
              int *d_evenArray,
              int *d_oddArray,
              spin *d_bEvenArray,
              spin *d_bOddArray,
              int *d_nn,
              int *d_accept,
              spin *d_phi_intermediate,
              double *d_aloc_comp,
              double *d_aloc_calc,
              double *d_rnd,
              double delta
              );
extern __global__ void alocal_gpu(spin *d_phi, spin *d_bArray, double *d_aloc);
/* extern __global__ spin magnet(void); */
/* extern __global__ double alocal(int idx); */
/* extern __global__ double alocal2(int idx, spin b); */
extern __global__ void compute_b(spin *d_phi, int *d_nn, spin *d_bArray);
extern __global__ void modify_spin(spin *d_phi, double *d_rnd, double delta, int ktrial);


#endif
