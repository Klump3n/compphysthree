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
     /* printf("rnd nums %f, %f, %f\n", rnd[0], rnd[1], rnd[2]); */
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

double metro_sweep_alt(double delta, int *evenArray, int *oddArray, double *rnd)
{
  int idx, acc;
  /* double *rnd; */

  int halfArrayLength = ceil((float) (nvol) / 2.);

  /* rnd=randgpu(nvol*3*NTRIAL); */

  acc=0;
  delta=2.0*delta;

  int evenIndex = 0;
  for (idx=0; idx<halfArrayLength; idx++)
    {
      evenIndex = evenArray[idx];
      acc+=update_spin(evenIndex,rnd,delta,NTRIAL);
      rnd+=(3*NTRIAL);
    }

  int oddIndex = 0;
  for (idx=0; idx<halfArrayLength; idx++)
    {
      oddIndex = oddArray[idx];
      acc+=update_spin(oddIndex,rnd,delta,NTRIAL);
      rnd+=(3*NTRIAL);
    }

  return ((double)(acc) / (double)(nvol*NTRIAL));
}
