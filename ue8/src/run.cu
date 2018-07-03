#define DEFINE_GLOBAL

#include "stat5.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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

   double *mm_array = (double *) malloc(nsweep*sizeof(double));
   double *s_array = (double *) malloc(nsweep*sizeof(double));

   bool thermalized = false;
   int thermCounter = 0;
   double varByExpMM = 0.0;
   int startIndex = 0;
   int N = 0;
   double factor_MM = 0.0,
     factor_S = 0.0;

   s=action();

   m=magnet();
   mm=cuCabs(m)*cuCabs(m);

   clear5(3, 2*nsweep);
   acc=0.0;
   iStart=seconds();
   for (i=1; i<=nsweep; i++)
   {
      acc=metro_sweep(delta);

      if (!(thermalized))
        {
          delta=tune_delta(acc,delta);
        }

      s=action();

      m=magnet();
      mm=cuCabs(m)*cuCabs(m);

      accum5(1, s);
      accum5(2, mm);
      accum5(3, aver5(2));

      varByExpMM = var5(2)/aver5(2)/aver5(2);
      if ((varByExpMM < .02) && i > 9)
        {
          thermCounter++;
        }

      /* when the varByExpMM was small 10 times */
      if ((thermCounter > 9) && !(thermalized))
        {
          printf("Thermalization took %d steps\n", i);
          printf("DELTA = %f\n", delta);
          printf("UPD\t S     \t\t |M|^2  \t avg(|M|^2)  \t std(|M|^2) \t tau(|M|^2) \t tau(S) \t DEF \n");
          thermalized = true;
          startIndex = i;
          i = 0;
        }

      if (thermalized)
        {
          N = i + startIndex;
          factor_MM = 2 * tauint5(2) * var5(2) / N / var5(3);
          /* mm_array[i] = mm; */
          /* s_array[i] = s; */
          /* printf("%d\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",i,acc,delta,s,cuCreal(m),cuCimag(m),mm, aver5(2), sigma5(2), tauint5(2), tauint5(1)); */
          printf("%d\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",i,s,mm, aver5(2), sigma5(2), tau5(2), tau5(1), factor_MM);
        }

   }

   printf("\n\n");
   printf("%d updates took %f sec.\n\n",nsweep,seconds()-iStart);

   free(lsize);
   free(nn[0]);
   free(nn);
   free(phi);
}
