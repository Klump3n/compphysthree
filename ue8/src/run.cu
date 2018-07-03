#define DEFINE_GLOBAL

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "stat5.h"
#include "global.h"
#include "randgpu.h"
#include "geom_pbc.h"
#include "metropolis.h"
#include "spin.h"
#include "common.h"


#define MIN_NARG 6

void usage(void)
{
   printf("Usage:\n\n");
   printf("  run <lambda> <kappa> <h> <phi0> <nsweep> <lsize1> [<lsize2> ...]\n\n");
   exit(0);
}

void init_phi(double phi0)
{
   int idx;

   if (phi0==0.0)
   {
      random_cnfg();
   }
   else
   {
      for (idx=0; idx<nvol; idx++)
      {
         phi[idx]=make_spin(phi0,0.0);
      }
   }
}

double tune_delta(double acc, double delta)
{
   if (acc<0.35)
      delta*=0.95;
   if (acc>0.45)
      delta*=1.05;

   return delta;
}

int main(int argc, char **argv)
{
   printf("%s Starting...\n", argv[0]);

   int i, nsweep;
   double acc, delta, s, iStart, mm, phi0, reh;
   spin m;

   if (argc<MIN_NARG+1)
      usage();

   // read parameters from command line
   lambda=atof(argv[1]);
   kappa=atof(argv[2]);
   reh=atof(argv[3]);
   phi0=atof(argv[4]);
   nsweep=atoi(argv[5]);

   ndim=argc-MIN_NARG;
   lsize=(int*)malloc((ndim+1)*sizeof(int));
   for (i=1; i<=ndim; i++)
   {
      lsize[i]=atoi(argv[i+MIN_NARG-1]);
   }

   delta=0.2;
   h=make_spin(reh,0.0);

   // print out parameters
   printf("Gittergroesse: %d",lsize[1]);
   for (i=2; i<=ndim; i++)
   {
      printf(" x %d",lsize[i]);
   }
   printf("\n\n");
   printf("nsweep = %d\n",nsweep);
   printf("lambda = %f\n",lambda);
   printf("kappa  = %f\n",kappa);
   printf("delta  = %f\n",delta);
   printf("h      = %f + I %f\n",cuCreal(h),cuCimag(h));
   if (phi0==0.0)
      printf("phi    = random\n");
   else
      printf("phi    = %f + I %f\n",phi0,0.0);
   printf("\n\n");

   // set up geometry
   geom_pbc();

   // allocate spins, set random values
   printf("Initalize spins... ");
   iStart=seconds();
   phi=(spin*)malloc(nvol*sizeof(spin));
   init_phi(phi0);
   printf("%f sec.\n\n",seconds()-iStart);

   bool thermalized = false;

   s=action();

   m=magnet();
   mm=cuCabs(m)*cuCabs(m);

   printf("UPD\t A       \t DELTA   \t S       \t RE(M)   \t IM(M)   \t |M|^2   \n");
   printf("%d\t %f\t %f\t %f\t %f\t %f\t %f\n",0,0.0,delta,s,cuCreal(m),cuCimag(m),mm);

   clear5(1, 500);
   acc=0.0;
   iStart=seconds();
   for (i=1; i<=nsweep; i++)
   {
      acc=metro_sweep(delta);
      delta=tune_delta(acc,delta);
      s=action();

      m=magnet();
      mm=cuCabs(m)*cuCabs(m);

      printf("%d\t %f\t %f\t %f\t %f\t %f\t %f\n",i,acc,delta,s,cuCreal(m),cuCimag(m),mm);

      /* accum5(1, s); */
      /* accum5(2, mm); */

   }

   printf("\n\n");
   printf("%d updates took %f sec.\n\n",nsweep,seconds()-iStart);

   free(lsize);
   free(nn[0]);
   free(nn);
   free(phi);
}
