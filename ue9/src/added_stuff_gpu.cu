#include "spin.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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

double gpu_sweep(
                 spin *d_phi,
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

  delta = 2.0 * delta;          /* consistency */

  int ntrial = 10;
  int acc_sum = 0;
  dim3 blockSize;

  int halfArrayLength = ceil((float) (nvol) / 2.);

  int *h_accept = (int *) calloc(halfArrayLength, sizeof(int));

  if (nvol/2 < 128)
    {
      blockSize.x = halfArrayLength;
    }
  else
    {
      blockSize.x = 128;
    }
  blockSize.y = 1;
  blockSize.z = 1;

  dim3 gridSize;
  gridSize.x = ceil((float) nvol / (float) 256);
  gridSize.y = 1;
  gridSize.z = 1;

  /* reset values of helper arrays */
  CHECK(cudaMemset(d_accept, 0, halfArrayLength*sizeof(int)));
  CHECK(cudaMemset(d_phi_intermediate, 0, nvol*sizeof(spin)));
  CHECK(cudaMemset(d_aloc_comp, 0.0, halfArrayLength*sizeof(double)));
  CHECK(cudaMemset(d_aloc_calc, 0.0, halfArrayLength*sizeof(double)));

  /* even sweep */
  /* b even */
  /**/
  compute_b_gpu<<<blockSize, gridSize>>>(d_phi, d_nn, d_bEvenArray);
  CHECK(cudaDeviceSynchronize());

  alocal_gpu<<<gridSize, blockSize>>>(d_phi, d_bEvenArray, d_aloc_comp);
  CHECK(cudaDeviceSynchronize());

  for (int i=0; i<ntrial; i++)
    {
      CHECK(cudaMemcpy(d_phi_intermediate, d_phi, nvol*sizeof(spin), cudaMemcpyDeviceToDevice));
      modify_spin<<<gridSize, blockSize>>>(d_phi, d_rnd, delta, i);
      CHECK(cudaDeviceSynchronize());
      alocal_gpu<<<gridSize, blockSize>>>(d_phi, d_bEvenArray, d_aloc_calc);
      CHECK(cudaDeviceSynchronize());
      gpu_comp_action <<<gridSize, blockSize>>> (d_phi, d_phi_intermediate, d_aloc_comp, d_aloc_calc, d_rnd, i, d_accept);
      CHECK(cudaDeviceSynchronize());
    }

  /* get d_accept and sum it all up */
  CHECK(cudaMemcpy(h_accept, d_accept, ((int)(nvol/2))*sizeof(int), cudaMemcpyDeviceToHost));
  for (int i=0; i<((int)(nvol/2)); i++)
    {
      acc_sum += h_accept[i];
    }

  /* end even */

  /* reset values of helper arrays */
  CHECK(cudaMemset(d_accept, 0, halfArrayLength*sizeof(int)));
  CHECK(cudaMemset(d_phi_intermediate, 0, nvol*sizeof(spin)));
  CHECK(cudaMemset(d_aloc_comp, 0.0, halfArrayLength*sizeof(double)));
  CHECK(cudaMemset(d_aloc_calc, 0.0, halfArrayLength*sizeof(double)));

  /* odd sweep */
  /* b odd */
  /**/
  compute_b_gpu<<<blockSize, gridSize>>>(d_phi, d_nn, d_bOddArray);
  CHECK(cudaDeviceSynchronize());

  alocal_gpu<<<gridSize, blockSize>>>(d_phi, d_bOddArray, d_aloc_comp);
  CHECK(cudaDeviceSynchronize());

  for (int i=0; i<ntrial; i++)
    {
      /* backup_spin <<<gridSize, blockSize>>> (d_phi, d_phi_intermediate); */
      /* CHECK(cudaDeviceSynchronize()); */
      CHECK(cudaMemcpy(d_phi_intermediate, d_phi, nvol*sizeof(spin), cudaMemcpyDeviceToDevice)); /* faster */
      modify_spin<<<gridSize, blockSize>>>(d_phi, d_rnd, delta, i);
      CHECK(cudaDeviceSynchronize());
      alocal_gpu<<<gridSize, blockSize>>>(d_phi, d_bOddArray, d_aloc_calc);
      CHECK(cudaDeviceSynchronize());
      gpu_comp_action <<<gridSize, blockSize>>> (d_phi, d_phi_intermediate, d_aloc_comp, d_aloc_calc, d_rnd, i, d_accept);
      CHECK(cudaDeviceSynchronize());
    }

  /* get d_accept and sum it all up */
  CHECK(cudaMemcpy(h_accept, d_accept, halfArrayLength*sizeof(int), cudaMemcpyDeviceToHost));
  for (int i=0; i<((int)(nvol/2)); i++)
    {
      acc_sum += h_accept[i];
    }

  /* end even */
  double percentage = (float) (acc_sum)/ (float) (ntrial*nvol);
  /* printf("%d of %d: %f\n", acc_sum, ntrial*nvol, percentage); */
  return percentage;

}

__global__
void backup_spin(spin *d_phi, spin *d_phi_intermediate)
{
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  d_phi_intermediate[idx] = d_phi[idx];
}

__global__
void compute_b_gpu(spin *d_phi, int *d_nn, spin *d_bArray)
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

__global__
void gpu_comp_action(spin *d_phi, spin *d_phi_intermediate, double *d_aloc_comp, double *d_aloc_calc, double *d_rnd, int ktrial, int *d_accept)
{
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

  /* NTRIAL = 10!!! */
  int rnd_idx = ktrial*10 + idx + 2;

  if (d_rnd[rnd_idx] < exp(-d_aloc_calc[idx] + d_aloc_comp[idx]))
    {
      d_accept[idx]++;
      d_aloc_comp[idx] = d_aloc_calc[idx];
    }
  else
    {
      d_phi[idx] = d_phi_intermediate[idx];
    }

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
