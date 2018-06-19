
#define DEFINE_GLOBAL

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include "global.h"

int main(void)
{
   int k;
   double r;
   cuDoubleComplex z = make_cuDoubleComplex(1.5,2.0);

   printf("  z    = %f + %f * i\n", cuCreal(z), cuCimag(z));
   r=cuCabs(z);
   printf(" |z|^2 = %.6f\n", r*r);
   z=cuCmul(z,cuConj(z));
   printf("  zz^* = %.6f + %.6f\n",  cuCreal(z), cuCimag(z));

   phi=(cuDoubleComplex*)malloc(10*sizeof(cuDoubleComplex));

   for (k=0; k<10; k++)
   {
      phi[k]=make_cuDoubleComplex((double)k,(double)(k*k));
      printf("%f + %f * i\n", cuCreal(phi[k]), cuCimag(phi[k]));
   }
}
