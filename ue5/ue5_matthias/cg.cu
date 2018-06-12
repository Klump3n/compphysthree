/*********************************************************************
cg.cu

Conjugate Gradient

double cg(double *x, double *r, int maxiter, double rel, int *status)
   Loest lineares Gleichungssystem
                  A*x = r
   maxiter: Maximale Anzahl der Iterationen
   rel:     Relative Reduktion der Residuumnorm
   status:  Wird bei Erfolg auf Anzahl der Iterationen gesetzt. Sonst <=0.

   Rueckgabewert ist die erreichte Residuumnorm.

**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "geometry.h"
#include "linalg.h"
#include "common.h"


#define DEBUG

double cg(double *x, double *r, int maxiter, double rel, int *status)
{
   int k;
   double ar,as,alpha,beta,rn,rnold,rn0;
   double *p,*s;

   s=(double*)malloc(npts*sizeof(double));
   p=(double*)malloc(npts*sizeof(double));

   memset(x,0,npts*sizeof(double));
   memset(s,0,npts*sizeof(double));

   rn=norm_sqr(r);
   rn0=rn;
   status[0]=0;
#ifdef DEBUG
   printf("Residuumnorm am Anfang: %e\n",sqrt(rn0));
#endif
   if (rn==0.0)
      return rn;

   assign_v2v(p,r);
   rel*=rel;
   k=0;

   while (k<maxiter)
   {
      laplace_2d(s,p);
      ar=vector_prod(p,r);
      as=vector_prod(p,s);
      alpha=ar/as;
      mul_add(x,alpha,p);
      mul_add(r,-alpha,s);
      rnold=rn;
      rn=norm_sqr(r);
      k+=1;
#ifdef DEBUG
      if (k % 10 == 0)
      {
         printf("Iter %d, rel. Residuumnorm: %e\n",k,sqrt(rn/rn0));
      }
#endif
      if ((rn/rn0)<rel)
      {
         break;
      }
      beta=rn/rnold;
      update_p(r,beta,p);
   }

#ifdef DEBUG
   printf("Rel. Residuumnorm nach %d Iterationen: %e\n",k,sqrt(rn/rn0));
#endif

   if ((rn/rn0<=rel) && (k<=maxiter))
      *status=k;
   if (rn/rn0>rel)
      *status=-1;

   free(s);
   free(p);

   return sqrt(rn);
}

double cg_gpu(double *x, double *r, int maxiter, double rel, int *status) {

  int k;
  double ar,as,alpha,beta,rn,rnold,rn0;
  double *d_p, *d_s, *d_r, *d_rr;

  CHECK(cudaMalloc((void **)&d_p, npts*sizeof(double)));
  CHECK(cudaMalloc((void **)&d_s, npts*sizeof(double)));
  CHECK(cudaMalloc((void **)&d_r, npts*sizeof(double)));
  CHECK(cudaMalloc((void **)&d_rr, npts*sizeof(double)));
  /* s=(double*)malloc(npts*sizeof(double)); */
  /* p=(double*)malloc(npts*sizeof(double)); */

  CHECK(cudaMemset(d_p, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_s, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_r, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_rr, 0, npts*sizeof(double)));
  /* memset(x,0,npts*sizeof(double)); */
  /* memset(s,0,npts*sizeof(double)); */

  /* r auf gpu laden */
  CHECK(cudaMemcpy(d_r, r, npts*sizeof(double), cudaMemcpyHostToDevice));

  rn = norm_sqr_gpu(r, d_rr, Nx, Ny);

  /* /\* unroll EC *\/ */
  /* unsigned int blockDimX = 128; /\* 32*k *\/ */
  /* unsigned int Nunroll = 8; */
  /* unsigned int gridDimX = Nx * Ny / Nunroll / blockDimX; */


  /* dim3 block2 (256,1); */
  /* int nblk = (npts + (block2.x*Nunroll) - 1)/(block2.x*Nunroll); */
  /* dim3 grid2 (nblk,1); */
  /* /\* rn=norm_sqr(r); *\/ */
  /* assign_v2v_gpu(d_rr, d_r); */
  /* vector_prod_pre_gpu<<<grid2, block2>>>(d_rr, d_r, Nx, Ny); */
  /* reduceUnrolling(d_rr, d_rr, npts); */
  /* rn = vector_add(d_rr, grid2.x); */
  /* /\* rn=norm_sqr(rr); *\/ */

  rn0=rn;
  status[0]=0;
#ifdef DEBUG
  printf("Residuumnorm am Anfang: %e\n",sqrt(rn0));
#endif
  if (rn==0.0)
    return rn;

  assign_v2v_gpu(d_p,d_r, Nx, Ny);
  rel*=rel;
  k=0;

  while (k<maxiter)
    {
      laplace_2d_gpu(d_s,d_p, Nx, Ny);
      /* ar=vector_prod(p,r); */
      /* as=vector_prod(p,s); */
      alpha=ar/as;
      mul_add_gpu(d_x,alpha,d_p, Nx, Ny);
      mul_add_gpu(d_r,-alpha,d_s, Nx, Ny);
      rnold=rn;
      rn=norm_sqr(r);
      k+=1;
#ifdef DEBUG
      if (k % 10 == 0)
        {
          printf("Iter %d, rel. Residuumnorm: %e\n",k,sqrt(rn/rn0));
        }
#endif
      if ((rn/rn0)<rel)
        {
          break;
        }
      beta=rn/rnold;
      update_p_gpu(r,beta,p, Nx, Ny);
    }

#ifdef DEBUG
  printf("Rel. Residuumnorm nach %d Iterationen: %e\n",k,sqrt(rn/rn0));
#endif

  if ((rn/rn0<=rel) && (k<=maxiter))
    *status=k;
  if (rn/rn0>rel)
    *status=-1;

  free(s);
  free(p);

  return sqrt(rn);
}
