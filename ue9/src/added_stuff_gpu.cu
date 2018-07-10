#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "spin.h"
#include "global.h"
#include "randgpu.h"
#include "added_stuff.h"
#include "added_stuff_gpu.h"

#define CHECK(call)                                             \
  {                                                             \
    const cudaError_t error = call;                             \
    if (error != cudaSuccess)                                   \
      {                                                         \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);  \
        fprintf(stderr, "code: %d, reason: %s\n", error,        \
                cudaGetErrorString(error));                     \
        exit(1);                                                \
      }                                                         \
  }

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
              )
{

  int ntrial = 10;
  int acc_sum = 0;
  int blockSize;

  if (nvol/2 < 128)
    {
      blockSize = (int) (nvol/2);
    }
  else
    {
      blockSize = 128;
    }

  int gridSize = ceil(nvol/256);

  /* even sweep */
  /* b even */
  /**/
  compute_b<<<gridSize, blockSize>>>(d_phi, d_nn, d_bEvenArray);
  alocal_gpu<<<gridSize, blockSize>>>(d_phi, d_bEvenArray, d_aloc_comp);

  for (int i=0; i<ntrial; i++)
    {
      CHECK(cudaMemcpy(d_phi_intermediate, d_phi, nvol*sizeof(spin), cudaMemcpyDeviceToDevice));
      modify_spin<<<gridSize, blockSize>>>(d_phi, d_rnd, delta, i);
      alocal_gpu<<<gridSize, blockSize>>>(d_phi, d_bEvenArray, d_aloc_calc);
      /* Vergleichsfunktion schreiben und das solls jewesen sein */
    }

  return acc_sum;
}

__global__
void modify_spin(spin *d_phi, double *d_rnd, double delta, int ktrial)
{
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

  /* NTRIAL = 10!!! */
  int rnd_idx = ktrial*10 + idx;

  d_phi[idx]=make_spin(
                       cuCreal(d_phi[idx])+delta*(d_rnd[rnd_idx + 0]-0.5),
                       cuCimag(d_phi[idx])+delta*(d_rnd[rnd_idx + 1]-0.5)
                       );

}

/* void gpu_update_spin() */
/* { */
/*   int k, acc; */
/*   spin tmp; */
/*   double a,ap; */

/*   acc=0; */
/*   a=alocal(idx); */
/*   for (k=0; k<ntrial; k++) */
/*     { */
/*       tmp=phi[idx]; */
/*       phi[idx]=make_spin(cuCreal(phi[idx])+del*(rnd[0]-0.5),cuCimag(phi[idx])+del*(rnd[1]-0.5)); */
/*       ap=alocal2(idx); */
/*       if (rnd[2]<exp(-ap+a)) */
/*         { */
/*           acc+=1; */
/*           a=ap; */
/*         } */
/*       else */
/*         phi[idx]=tmp; */

/*       rnd+=3; */
/*     } */

/*   return acc; */

/* } */


/* __global__ */
/* spin magnet(void) */
/* { */
/*   int idx; */
/*   spin m; */

/*   m=make_spin(0.0,0.0); */
/*   for (idx=0; idx<nvol; idx++) */
/*     { */
/*       m=cuCadd(m,phi[idx]); */
/*     } */
/*   m=make_spin(cuCreal(m)/nvol,cuCimag(m)/nvol); */

/*   return m; */
/* } */

__global__
void compute_b(spin *d_phi, int *d_nn, spin *d_bArray)
{
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

  spin tmpc;

  d_bArray[idx] = devH;
  for (int k=1; k<=devNdim; k++)
    {
      tmpc = cuCadd(d_phi[d_nn[k*devNvol+ idx]], d_phi[d_nn[(devNdim+k)*devNvol + idx]]);
      d_bArray[idx] = make_cuDoubleComplex(cuCreal(d_bArray[idx])+devKappa*cuCreal(tmpc), cuCimag(d_bArray[idx])+devKappa*cuCimag(tmpc));
    }
}

// use pre-computed b
__global__
void alocal_gpu(spin *d_phi, spin *d_bArray, double *d_aloc)
{
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

  double a,tmp;

  tmp=cuCreal(d_phi[idx])*cuCreal(d_phi[idx]) + cuCimag(d_phi[idx])*cuCimag(d_phi[idx]);
  a=2.0*(cuCreal(d_bArray[idx])*cuCreal(d_phi[idx])+cuCimag(d_bArray[idx])*cuCimag(d_phi[idx]))-tmp;
  tmp-=1.0;
  a-=devLambda*tmp*tmp;

  d_aloc[idx] = -a;
}

/* // compute b */
/* __global__ */
/* double alocal() */
/* { */
/*   unsigned int tid = threadIdx.x; */
/*   unsigned int gridSize = blockDim.x*gridDim.x; */
/*   unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x; */

/*   spin b = compute_b(idx); */

/*   return alocal2(idx, b); */
/* } */
