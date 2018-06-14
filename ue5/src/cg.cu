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
      /* printf("cpu norm laplace = %f\n", norm_sqr(s)); */
      ar=vector_prod(p,r);
      /* printf("cpu ar = %f\n", ar); */
      as=vector_prod(p,s);
      /* printf("cpu as = %f\n", as); */
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
  double *d_p, *d_s, *d_r, *d_x, *d_intmed1, *d_intmed2, *d_intmed3;

  CHECK(cudaMalloc((void **)&d_p, npts*sizeof(double)));
  CHECK(cudaMalloc((void **)&d_s, npts*sizeof(double)));
  CHECK(cudaMalloc((void **)&d_r, npts*sizeof(double)));
  CHECK(cudaMalloc((void **)&d_x, npts*sizeof(double)));

  /* allocate generous space */
  CHECK(cudaMalloc((void **)&d_intmed1, npts*sizeof(double)));
  CHECK(cudaMalloc((void **)&d_intmed2, npts*sizeof(double)));
  CHECK(cudaMalloc((void **)&d_intmed3, npts*sizeof(double)));

  CHECK(cudaMemset(d_p, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_s, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_r, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_x, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_intmed1, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_intmed2, 0, npts*sizeof(double)));
  CHECK(cudaMemset(d_intmed3, 0, npts*sizeof(double)));

  double iStart = seconds();

  /* r auf gpu laden */
  CHECK(cudaMemcpy(d_r, r, npts*sizeof(double), cudaMemcpyHostToDevice));

  rn = norm_sqr(r);             /* this is quicker than the gpu version */
  rn0=rn;
  status[0]=0;
#ifdef DEBUG
  printf("Residuumnorm am Anfang: %e\n",sqrt(rn0));
#endif
  if (rn==0.0)
    return rn;

  assign_v2v_gpu<<<grid, block>>>(d_p,d_r, Nx, Ny);
  /* CHECK(cudaDeviceSynchronize()); */

  rel*=rel;
  k=0;

  while (k<maxiter)
    {
      laplace_2d_gpu<<<grid, block>>>(d_s,d_p, Nx, Ny);
      /* CHECK(cudaDeviceSynchronize()); */

      vector_prod_gpu(&ar, d_p, d_r, d_intmed1, d_intmed2, d_intmed3, Nx, Ny);
      vector_prod_gpu(&as, d_p, d_s, d_intmed1, d_intmed2, d_intmed3, Nx, Ny);

      alpha=ar/as;
      mul_add_gpu<<<grid, block>>>(d_x,alpha,d_p, Nx, Ny);
      /* CHECK(cudaDeviceSynchronize()); */
      mul_add_gpu<<<grid, block>>>(d_r,-alpha,d_s, Nx, Ny);
      /* CHECK(cudaDeviceSynchronize()); */
      rnold=rn;
      norm_sqr_gpu(&rn, d_r, d_intmed1, d_intmed2, Nx, Ny);

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
      update_p_gpu<<<grid, block>>>(d_r,beta,d_p, Nx, Ny);
      /* CHECK(cudaDeviceSynchronize()); */
    }

#ifdef DEBUG
  printf("Rel. Residuumnorm nach %d Iterationen: %e\n",k,sqrt(rn/rn0));
#endif

  if ((rn/rn0<=rel) && (k<=maxiter))
    *status=k;
  if (rn/rn0>rel)
    *status=-1;

  return sqrt(rn);





/*   assign_v2v(p,r); */
/*   rel*=rel; */
/*   k=0; */

/*   while (k<maxiter) */
/*     { */
/*       laplace_2d(s,p); */
/*       ar=vector_prod(p,r); */
/*       as=vector_prod(p,s); */
/*       alpha=ar/as; */
/*       mul_add(x,alpha,p); */
/*       mul_add(r,-alpha,s); */
/*       rnold=rn; */
/*       rn=norm_sqr(r); */
/*       k+=1; */
/* #ifdef DEBUG */
/*       if (k % 10 == 0) */
/*         { */
/*           printf("Iter %d, rel. Residuumnorm: %e\n",k,sqrt(rn/rn0)); */
/*         } */
/* #endif */
/*       if ((rn/rn0)<rel) */
/*         { */
/*           break; */
/*         } */
/*       beta=rn/rnold; */
/*       update_p(r,beta,p); */
/*     } */

/* #ifdef DEBUG */
/*   printf("Rel. Residuumnorm nach %d Iterationen: %e\n",k,sqrt(rn/rn0)); */
/* #endif */

/*   if ((rn/rn0<=rel) && (k<=maxiter)) */
/*     *status=k; */
/*   if (rn/rn0>rel) */
/*     *status=-1; */

/*   free(s); */
/*   free(p); */

/*   return sqrt(rn); */




}
