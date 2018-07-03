#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "spin.h"
#include "randgpu.h"

#define NTRIAL 10

int update_spin(int idx, double* rnd, double del, int ntrial)
{
   int k, acc;
   spin tmp;
   double a,ap;

   acc=0;
   a=alocal(idx);
   for (k=0; k<ntrial; k++)
   {
      tmp=phi[idx];
      phi[idx]=make_spin(cuCreal(phi[idx])+del*(rnd[0]-0.5),cuCimag(phi[idx])+del*(rnd[1]-0.5));
      ap=alocal2(idx);
      if (rnd[2]<exp(-ap+a))
      {
         acc+=1;
         a=ap;
      }
      else
         phi[idx]=tmp;

      rnd+=3;
   }

   return acc;
}

double metro_sweep(double delta)
{
   int idx, acc;
   double *rnd;

   rnd=randgpu(nvol*3*NTRIAL);

   acc=0;
   delta=2.0*delta;
   for (idx=0; idx<nvol; idx++)
   {
      acc+=update_spin(idx,rnd,delta,NTRIAL);
      rnd+=(3*NTRIAL);
   }

   return ((double)(acc) / (double)(nvol*NTRIAL));
}
