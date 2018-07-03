#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "spin.h"
#include "global.h"
#include "randgpu.h"

static spin b;

double action(void)
{
   int idx,k;
   double act,tmp,ka2;
   spin h2;

#ifdef DEBUG
   printf("lambda: %f, kappa: %f, h: %f + I %f\n",lambda,kappa,cuCreal(h),cuCimag(h));
#endif

   act=0.0;
   ka2=2.0*kappa;
   h2=make_spin(cuCreal(h)*2.0,cuCimag(h)*2.0);
   for (idx=0; idx<nvol; idx++)
   {
      tmp=cuCreal(phi[idx])*cuCreal(phi[idx]) + cuCimag(phi[idx])*cuCimag(phi[idx]);
      act+=tmp;
      tmp-=1.0;
      act+=lambda*tmp*tmp;
      tmp=0.0;
      for (k=1; k<=ndim; k++)
      {
         tmp+=(cuCreal(phi[idx])*cuCreal(phi[nn[k][idx]]) + cuCimag(phi[idx])*cuCimag(phi[nn[k][idx]]));
      }
      act-=ka2*tmp;
      act-=(cuCreal(phi[idx])*cuCreal(h2) + cuCimag(phi[idx])*cuCimag(h2));
   }

   return act;
}

void random_cnfg(void)
{
   int idx;
   double *rnd;

   rnd=randgpu(2*nvol);

   for (idx=0; idx<nvol; idx++)
   {
      phi[idx]=make_spin(rnd[2*idx]-0.5,(rnd[2*idx+1]-0.5));
   }
}

spin magnet(void)
{
   int idx;
   spin m;

   m=make_spin(0.0,0.0);
   for (idx=0; idx<nvol; idx++)
   {
      m=cuCadd(m,phi[idx]);
   }
   m=make_spin(cuCreal(m)/nvol,cuCimag(m)/nvol);

   return m;
}

void compute_b(int idx)
{
   int k;
   spin tmpc;

   b=h;
   for (k=1; k<=ndim; k++)
   {
      tmpc=cuCadd(phi[nn[k][idx]],phi[nn[ndim+k][idx]]);
      b=make_cuDoubleComplex(cuCreal(b)+kappa*cuCreal(tmpc),cuCimag(b)+kappa*cuCimag(tmpc));
   }
}

// use pre-computed b
double alocal2(int idx)
{
   double a,tmp;

   tmp=cuCreal(phi[idx])*cuCreal(phi[idx]) + cuCimag(phi[idx])*cuCimag(phi[idx]);
   a=2.0*(cuCreal(b)*cuCreal(phi[idx])+cuCimag(b)*cuCimag(phi[idx]))-tmp;
   tmp-=1.0;
   a-=lambda*tmp*tmp;

   #ifdef DEBUG
      printf("b: %f + %f i, a: %f\n",cuCreal(b),cuCimag(b),a);
   #endif

   return -a;
}

// compute b
double alocal(int idx)
{
   compute_b(idx);

   return alocal2(idx);
}
