#ifndef ADDED_STUFF_GPU_H
#define ADDED_STUFF_GPU_H

#include <cuComplex.h>

typedef cuDoubleComplex spin;

#define make_spin(a,b) (                          \
                        make_cuDoubleComplex(a,b) \
                         )

double gpu_sweep(spin *d_phi,
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
double cpu_sweep(
                 spin *h_phi,
                 int *h_evenArray,
                 int *h_oddArray,
                 spin *h_bEvenArray,
                 spin *h_bOddArray,
                 int *h_nn,
                 int *h_accept,
                 spin *h_phi_intermediate,
                 double *h_aloc_comp,
                 double *h_aloc_calc,
                 double *h_rnd,
                 double delta
                 );
__global__ void alocal_gpu(spin *d_phi, spin *d_bArray, double *d_aloc);
/* extern __global__ spin magnet(void); */
/* extern __global__ double alocal(int idx); */
/* extern __global__ double alocal2(int idx, spin b); */
__global__ void compute_b_gpu(spin *d_phi, int *d_nn, spin *d_bArray);
__global__ void modify_spin(spin *d_phi, double *d_rnd, double delta, int ktrial);
__global__ void gpu_comp_action(spin *d_phi, spin *d_phi_intermediate, double *d_aloc_comp, double *d_aloc_calc, double *d_rnd, int ktrial, int *d_accept);
/* __global__ void backup_spin(spin *d_phi, spin *d_phi_intermediate); */
__global__ void reduceUnrolling (int *g_idata, int *g_odata, unsigned int n);

#endif
