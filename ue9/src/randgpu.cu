/*

   Zufallszahlen von der GPU
      Siehe: https://docs.nvidia.com/cuda/curand/index.html

   double* randgpu(unsigned int N)
      Erzeugt N gleichvertielte Zufallszahlen in (0,1] auf der GPU und
      kopiert diese in den Host-Speicher. Host-Speicher wir intern
      alloziert und gemanagt. Rückgabewert ist der Zeiger auf den
      Host-Speicher.

*/

#include <stdio.h>
#include <stdlib.h>
#include <curand_kernel.h> // CURAND Bibliothek header-Datei
#include "common.h"
#include "randgpu.h"

#define BLOCKSIZE 256

//#define DEBUG

static curandState *d_states=NULL;
static double *d_rnd=NULL;
static double *h_rnd=NULL;
static unsigned int gridsize=0;

__global__ void init_curand(unsigned int seed, curandState *states)
{
   unsigned int tid  = threadIdx.x + blockDim.x*blockIdx.x;

   // Initialisieren des Zufallszahlen-Generators
   // Der 'state' (Status) wird für jeden Thread unabhängig gespeichert
   curand_init(seed, tid, 0, &states[tid]);
}

__global__ void rnd_gpu(double *d_rnd, curandState *states)
{
   unsigned int tid  = threadIdx.x + blockDim.x*blockIdx.x;;

   d_rnd[tid]=curand_uniform_double(&states[tid]);
}

static void alloc_random(unsigned int grd)
{
   int nth;

   nth=grd*BLOCKSIZE;
   h_rnd=(double*)malloc(nth*sizeof(double));
   CHECK(cudaMalloc((void**)&d_rnd,nth*sizeof(double)));
   CHECK(cudaMalloc((void**)&d_states,nth*sizeof(curandState)));
#ifdef DEBUG
      printf("Initialize random number generator.\n");
#endif
   init_curand<<<grd,BLOCKSIZE>>>((unsigned int)seconds(),d_states);
   CHECK(cudaDeviceSynchronize());
   CHECK(cudaGetLastError());
   gridsize=grd;
}

double* randgpu(unsigned int N)
{
   unsigned int grd;

   grd=(N+BLOCKSIZE-1)/BLOCKSIZE;

   if (d_states==NULL)
      alloc_random(grd);

   if (gridsize<grd)
   {
      #ifdef DEBUG
         printf("reallocate buffers, max. grid: %d, block: %d, N: %d, grid: %d\n",gridsize,BLOCKSIZE,N,grd);
      #endif
      free(h_rnd);
      CHECK(cudaFree(d_rnd));
      CHECK(cudaFree(d_states));
      alloc_random(grd);
   }

#ifdef DEBUG
   printf("max. grid: %d, block: %d, N: %d, grid: %d\n",gridsize,BLOCKSIZE,N,grd);
#endif

   rnd_gpu<<<grd,BLOCKSIZE>>>(d_rnd,d_states);
   CHECK(cudaDeviceSynchronize());
   CHECK(cudaMemcpy(h_rnd, d_rnd, N*sizeof(double), cudaMemcpyDeviceToHost));
   CHECK(cudaGetLastError());

   return h_rnd;
}

double* randgpu_device_ptr(unsigned int N)
{
  unsigned int grd;

  grd=(N+BLOCKSIZE-1)/BLOCKSIZE;

  if (d_states==NULL)
    alloc_random(grd);

  if (gridsize<grd)
    {
#ifdef DEBUG
      printf("reallocate buffers, max. grid: %d, block: %d, N: %d, grid: %d\n",gridsize,BLOCKSIZE,N,grd);
#endif
      free(h_rnd);
      CHECK(cudaFree(d_rnd));
      CHECK(cudaFree(d_states));
      alloc_random(grd);
    }

#ifdef DEBUG
  printf("max. grid: %d, block: %d, N: %d, grid: %d\n",gridsize,BLOCKSIZE,N,grd);
#endif

  rnd_gpu<<<grd,BLOCKSIZE>>>(d_rnd,d_states);
  CHECK(cudaDeviceSynchronize());

  return d_rnd;
}
