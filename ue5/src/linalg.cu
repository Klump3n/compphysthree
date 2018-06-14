/*********************************************************************
linalg.cu

Funktionen die Vektor-Operationen und Matrix-Vektor-Opertionen implementieren:

void random_vector(double *p)
   Setzt p im Inneren auf zufaellige Werte

__global__ void reduceUnrolling (double *g_idata, double *g_odata, unsigned int n)
   Reduktion auf der GPU. Berechnet die Teilsumme auf jedem Block mit vorheriger
   maximaler Anzahl von 'seriellen' Additionen pro Thread.

double norm_sqr(double *v)
   Quadrat der Norm von v.

double vector_prod(double *v, double *w)
   Vektorprodukt v^T * w

void assign_v2v(double *v, double *w)
__global__ void assign_v2v_gpu(double *v, double *w, int nx, int ny)
   Zuweisung v = w

void mul_add(double *v, double a, double *w)
__global__ void mul_add_gpu(double *v, double a, double *w, int nx, int ny)
   Multiplikation von w mit Skalar a und Addition zu v. Ergebnis in v.
                  v = v + a*w

void update_p(double *r, double b, double *p)
__global__ void update_p_gpu(double *r, double b, double *p, int nx, int ny)
   Multiplikation von p mit Skalar b und Addition zu r. Ergebnis in p.
                  p = r + b*p

void laplace_2d(double *w, double *v)
__global__ void laplace_2d_gpu(double *w, double *v, int nx, int ny)
   Anwendung des 2-D Laplace-Operator A=-\Delta mit Dirichlet-Randbedingungen
                  w = A*v

**********************************************************************/
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "global.h"
#include "geometry.h"
#include "linalg.h"
#include "common.h"

/*
   Der Vektor p wird im Inneren auf zufaellige Werte gesetzt
*/
void random_vector(double *p)
{
   int idx;

   for(idx = 0; idx < npts; idx++)
   {
      if (active[idx])
         p[idx] = (double)(rand() & 0xFF ) / 10.0;
   }
}

__global__ void reduceUnrolling (double *g_idata, double *g_odata, unsigned int n)
{
    // set thread ID
    unsigned int tid = threadIdx.x;
    unsigned int gridSize = blockDim.x*2*gridDim.x;
    unsigned int idx = blockIdx.x * blockDim.x * 2 + threadIdx.x;

    // unroll as many as possible
    unsigned int nunroll=n/gridSize, k;
    unsigned int i=idx+nunroll*gridSize;
    double sum=0.0;
    if (i<n)
        sum += g_idata[i];
    if (i+blockDim.x<n)
        sum += g_idata[i+blockDim.x];
    for (k=1; k<=nunroll; k++)
    {
        i -= gridSize;
        sum += g_idata[i] + g_idata[i+blockDim.x];
    }
    g_idata[idx] = sum;

    __syncthreads();

    // in-place reduction in global memory
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2)
    {
        if (tid < stride)
        {
            g_idata[idx] += g_idata[idx + stride];
        }

        // synchronize within threadblock
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = g_idata[idx];
}

double norm_sqr(double *v)
{
   int idx;
   double r=0.0;
   for (idx=0; idx<npts; idx++)
   {
      r+=v[idx]*v[idx];
   }
   return r;
}
/* this is super slow */
void norm_sqr_gpu(double *res, double *d_r, double *d_intmed1, double *d_intmed2, int nx, int ny) {

  /* set intmeds to 0 */
  CHECK(cudaMemset(d_intmed1, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_intmed2, 0, npts*sizeof(double)));

  /* square all entries of the vector */
  vector_square_entries_gpu<<<grid, block>>>(d_intmed1, d_r, nx, ny);
  /* CHECK(cudaDeviceSynchronize()); */

  /* add as much as possible */
  reduceUnrolling<<<grid2, block2>>>(d_intmed1, d_intmed2, npts);
  /* CHECK(cudaDeviceSynchronize()); */

  /* add the rest too */
  reduceUnrolling<<<grid3, block3>>>(d_intmed2, d_intmed1, nblk);
  /* CHECK(cudaDeviceSynchronize()); */

  CHECK(cudaMemcpy(res, d_intmed1, 1*sizeof(double),cudaMemcpyDeviceToHost));

  /* CHECK(cudaDeviceSynchronize()); */
}

__global__ void alt_vect_prod(double *d_res, double *d_v, double *d_w, int nx, int ny) {
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x + 1;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y + 1;
  unsigned int idx = iy * (nx+2) + ix;

  if (ix<=nx && iy<=ny)
    {
      *d_res += d_v[idx]*d_w[idx];
    }
}
void vector_prod_gpu(double *res, double *d_v, double *d_w, double *d_intmed1, double *d_intmed2, double *d_intmed3, int nx, int ny) {

  /* set intmeds to 0 */
  CHECK(cudaMemset(d_intmed1, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_intmed2, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_intmed3, 0, npts*sizeof(double)));

  /* multiply all the vector entries and store the results in d_intmed1*/
  vector_multiply_entries_gpu<<<grid, block>>>(d_intmed1, d_v, d_w, nx, ny);
  /* CHECK(cudaDeviceSynchronize()); */

  /* sum up all the entries */
  /* add as much as possible */
  reduceUnrolling<<<grid2, block2>>>(d_intmed1, d_intmed2, npts);
  /* CHECK(cudaDeviceSynchronize()); */

  /* add the rest too */
  reduceUnrolling<<<grid3, block3>>>(d_intmed2, d_intmed3, nblk);
  /* CHECK(cudaDeviceSynchronize()); */

  CHECK(cudaMemcpy(res, d_intmed2, sizeof(double),cudaMemcpyDeviceToHost));

  /* CHECK(cudaDeviceSynchronize()); */
}

double vector_prod(double *v, double *w)
{
   int idx;
   double r=0.0;
   for (idx=0; idx<npts; idx++)
   {
      r+=v[idx]*w[idx];
   }
   return r;
}

double vector_add(double *v, const int n) {
  int idx;
  double r=0.0;
  for (idx=0; idx<n; idx++)
    {
      r+=v[idx];
    }
  return r;
}

__global__ void vector_square_entries_gpu(double *v, double *w, int nx, int ny) {
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x + 1;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y + 1;
  unsigned int idx = iy * (nx+2) + ix;

  if (ix<=nx && iy<=ny)
    {
      v[idx] = w[idx]*w[idx];
    }
}

__global__ void vector_multiply_entries_gpu(double *res, double *v, double *w, int nx, int ny)
{
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x + 1;
  unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y + 1;
  unsigned int idx = iy * (nx+2) + ix;

  if (ix<=nx && iy<=ny)
    {
      res[idx] = v[idx] * w[idx];
    }
}

void assign_v2v(double *v, double *w)
{
   int idx;
   for (idx=0; idx<npts; idx++)
   {
      v[idx]=w[idx];
   }
}

__global__ void assign_v2v_gpu(double *v, double *w, int nx, int ny)
{
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x + 1;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y + 1;
   unsigned int idx = iy * (nx+2) + ix;

   if (ix<=nx && iy<=ny)
   {
      v[idx]=w[idx];
   }
}

void mul_add(double *v, double a, double *w)
{
   int idx;
   for (idx=0; idx<npts; idx++)
   {
      v[idx]+=a*w[idx];
   }
}

__global__ void mul_add_gpu(double *v, double a, double *w, int nx, int ny)
{
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x + 1;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y + 1;
   unsigned int idx = iy * (nx+2) + ix;

   if (ix<=nx && iy<=ny)
   {
      v[idx]+=a*w[idx];
   }
}

void update_p(double *r, double b, double *p)
{
   int idx;
   for (idx=0; idx<npts; idx++)
   {
      p[idx]=r[idx]+b*p[idx];
   }
}

__global__ void update_p_gpu(double *r, double b, double *p, int nx, int ny)
{
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x + 1;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y + 1;
   unsigned int idx = iy * (nx+2) + ix;

   if (ix<=nx && iy<=ny)
   {
      p[idx]=r[idx]+b*p[idx];
   }
}

/*
   2D Laplace-Operator A=-\Delta, Dirichlet-Randbedingungen,
   multipliziert mit Vektor v:

                  w = A*v

   Wirkt nur auf innere/aktive Punkte. Aeussere Punkte bleiben unveraendert.
*/
void laplace_2d(double *w, double *v)
{
   int idx;
   for (idx=0; idx<npts; idx++)
   {
      if (active[idx])
         w[idx]=4.0*v[idx] - v[idx+1] - v[idx-1] - v[idx+Nx+2] - v[idx-Nx-2];
   }
}

__global__ void laplace_2d_gpu(double *w, double *v, int nx, int ny)
{
   unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x + 1;
   unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y + 1;
   unsigned int idx = iy * (nx+2) + ix;

   if (ix<=nx && iy<=ny)
   {
      w[idx]=4.0*v[idx] - v[idx+1] - v[idx-1] - v[idx+nx+2] - v[idx-nx-2];
   }
}
