#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "common.h"
#include "spin.h"
#include "randgpu.h"
#include "metropolis.h"

#define BLOCKSIZE 256
//#define DEBUG

/* Host and Device */

// compute b
__host__ __device__ spin blocal(spin *p, int idx, int *nni, int nd, int nv, double ka, spin hf)
{
   int k,in;
   double r,i;

   r=cuCreal(hf);
   i=cuCimag(hf);
   for (k=1; k<=nd; k++)
   {
      in=nni[k*nv+idx];
      r+=cuCreal(p[in]);
      i+=cuCimag(p[in]);
      in=nni[(k+nd)*nv+idx];
      r+=cuCreal(p[in]);
      i+=cuCimag(p[in]);
   }

   return make_spin(ka*r,ka*i);
}

__host__ __device__ double alocal3(spin* p, spin b, double la)
{
   double a,tmp;

   tmp=cuCreal(*p)*cuCreal(*p) + cuCimag(*p)*cuCimag(*p);
   a=2.0*(cuCreal(b)*cuCreal(*p)+cuCimag(b)*cuCimag(*p))-tmp;
   tmp-=1.0;
   a-=la*tmp*tmp;

   #ifdef DEBUG
      printf("b: %f + %f i, a: %f\n",cuCreal(b),cuCimag(b),a);
   #endif

   return -a;
}

__host__ __device__ int update_spin2(spin* p, spin b, double* rnd, double del, double la, int ntrial)
{
   int k, acc;
   spin tmp;
   double a,ap;

   acc=0;
   a=alocal3(p,b,la);
   for (k=0; k<ntrial; k++)
   {
      tmp=*p;
      *p=make_spin(cuCreal(*p)+del*(rnd[0]-0.5),cuCimag(*p)+del*(rnd[1]-0.5));
      ap=alocal3(p,b,la);
      if (rnd[2]<exp(-ap+a))
      {
         acc+=1;
         a=ap;
      }
      else
         *p=tmp;

      rnd+=3;
   }

   return acc;
}

/* end of Host and Device */

/* Device-Version */

int *acc, *d_eo, *d_acc, *d_nn;
spin *d_phi;
int gridsize, blocksize;
__device__ int *devEo, *devAcc;
__device__ double devLambda, devKappa;
__device__ int devNdim, devNvol, *devNn;
__device__ spin devH;
__device__ spin *devPhi;

__global__ void update_eo(int ieo, curandState *states, double del, int ntrial)
{
   unsigned int tid  = threadIdx.x + blockDim.x*blockIdx.x;

   int idx, k;
   double rnd[3*NTRIAL];
   spin b;

   if (tid>devNvol/2-1) return;

   idx=devEo[ieo*(devNvol/2)+tid];
   devAcc[idx]=0;

#ifdef DEBUG
   printf("idx: %d, tid: %d\n",idx,tid);
   printf("devNvol: %d, devNdim: %d, devKappa: %f, devLambda: %f\n",devNvol,devNdim,devKappa,devLambda);
#endif

   b=blocal(devPhi,idx,devNn,devNdim,devNvol,devKappa,devH);

   for (k=0; k<3*ntrial; k++)
      rnd[k]=curand_uniform_double(&states[tid]);

   devAcc[idx]=update_spin2(&devPhi[idx],b,rnd,del,devLambda,ntrial);
}

__global__ void update_eo_rnd(int ieo, double *rnd, double del, int ntrial)
{
   unsigned int tid  = threadIdx.x + blockDim.x*blockIdx.x;

   int idx;
   spin b;

   if (tid>devNvol/2-1) return;

   idx=devEo[ieo*(devNvol/2)+tid];
   devAcc[idx]=0;

#ifdef DEBUG
   printf("idx: %d, tid: %d\n",idx,tid);
   printf("devNvol: %d, devNdim: %d, devKappa: %f, devLambda: %f\n",devNvol,devNdim,devKappa,devLambda);
#endif

   // compute B
   b=blocal(devPhi,idx,devNn,devNdim,devNvol,devKappa,devH);

   rnd+=(idx*3*ntrial);
   devAcc[idx]=update_spin2(&devPhi[idx],b,rnd,del,devLambda,ntrial);
}

int* update_spin_eo_gpu(int ieo, double* rnd, double del, int ntrial)
{
   curandState *states;

#ifdef DEBUG
   printf("update_spin_gpu:\n ieo: %d, <<<%d, %d>>>\n",ieo,gridsize,blocksize);
#endif

   if (rnd==NULL)
   {
      states=get_curand_states(gridsize*blocksize);
      update_eo<<<gridsize,blocksize>>>(ieo, states, del, ntrial);
   }
   else
   {
      update_eo_rnd<<<gridsize,blocksize>>>(ieo, rnd, del, ntrial);
   }
   CHECK(cudaDeviceSynchronize());

   return d_acc;
}

void init_gpu(void)
{
   // Parameter
   CHECK(cudaMemcpyToSymbol(devLambda, &lambda, sizeof(double)));
   CHECK(cudaMemcpyToSymbol(devKappa, &kappa, sizeof(double)));
   CHECK(cudaMemcpyToSymbol(devNdim, &ndim, sizeof(int)));
   CHECK(cudaMemcpyToSymbol(devNvol, &nvol, sizeof(int)));
   CHECK(cudaMemcpyToSymbol(devH, &h, sizeof(spin)));

   // naechste Nachbarn
   CHECK(cudaMalloc((void**)&d_nn,nvol*(2*ndim+1)*sizeof(int)));
   CHECK(cudaMemcpy(d_nn, nn[0], nvol*(2*ndim+1)*sizeof(int), cudaMemcpyHostToDevice));
   CHECK(cudaMemcpyToSymbol(devNn, &d_nn, sizeof(int*)));

   // Even-odd Liste
   CHECK(cudaMalloc((void**)&d_eo,nvol*sizeof(int)));
   CHECK(cudaMemcpy(d_eo, eo[0], nvol*sizeof(int), cudaMemcpyHostToDevice));
   CHECK(cudaMemcpyToSymbol(devEo, &d_eo, sizeof(int*)));

   // Spin-Feld
   CHECK(cudaMalloc((void**)&d_phi,nvol*sizeof(spin)));
   CHECK(cudaMemcpyToSymbol(devPhi, &d_phi, sizeof(spin*)));
   CHECK(cudaMemcpy(d_phi, phi, nvol*sizeof(spin), cudaMemcpyHostToDevice));

   // Feld fuer Akzeptanz
   acc = (int *) malloc(nvol*sizeof(int));
   CHECK(cudaMalloc((void**)&d_acc,nvol*sizeof(int)));
   CHECK(cudaMemcpyToSymbol(devAcc, &d_acc, sizeof(int*)));

   // Execution-Konfiguration
   blocksize=BLOCKSIZE;
   gridsize=(nvol/2+blocksize-1)/blocksize;
}

void gather_phi(void)
{
   CHECK(cudaMemcpy(phi, d_phi, nvol*sizeof(spin), cudaMemcpyDeviceToHost));
}

double metro_sweep_gpu(double delta, double *rnd)
{
   int acca, idx;

   // even
   update_spin_eo_gpu(0,rnd,2.0*delta,NTRIAL);
   // odd
   update_spin_eo_gpu(1,rnd,2.0*delta,NTRIAL);

   // Akzeptanz
   CHECK(cudaMemcpy(acc, d_acc, nvol*sizeof(int), cudaMemcpyDeviceToHost));
   acca=0;
   for (idx=0; idx<nvol; idx++)
      acca+=acc[idx];

   // Feld kopieren
   gather_phi();

   return ((double)(acca) / (double)(nvol*NTRIAL));
}

/* end of Device-Version */

/* Host Version */

int update_spin(int idx, double* rnd, double del, int ntrial)
{
   spin b;

   b=compute_b(idx);
   return update_spin2(&phi[idx],b,rnd,del,lambda,ntrial);
}

double metro_sweep_cpu(double delta, double *rnd)
{
   int idx, acc, ix;

   acc=0;
   delta=2.0*delta;
   for (idx=0; idx<nvol/2; idx++)
   {
      ix=eo[0][idx];
      acc+=update_spin(ix,&rnd[ix*3*NTRIAL],delta,NTRIAL);
   }
   for (idx=0; idx<nvol/2; idx++)
   {
      ix=eo[1][idx];
      acc+=update_spin(ix,&rnd[ix*3*NTRIAL],delta,NTRIAL);
   }

   return ((double)(acc) / (double)(nvol*NTRIAL));
}

/* end of Host-Version */

double metro_sweep(double delta, double *rnd, int gpu)
{
   if (gpu>0)
      return metro_sweep_gpu(delta,rnd); // Run on Device
   else
      return metro_sweep_cpu(delta,rnd); // Run on Host
}
