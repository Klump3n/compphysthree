#define DEFINE_GLOBAL

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
#include "stat5.h"

#define MIN_NARG 7
#define NTHERM 1000

void usage(void)
{
   printf("Usage:\n\n");
   printf("  run <lambda> <kappa> <h> <phi0> <nsweep> <use_gpu> <lsize1> [<lsize2> ...]\n\n");
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

double* rand4sweep(int use_gpu)
{
   double *rnd;

   // Zufallszahlen
   if (use_gpu==1)
      rnd=devRandgpu(nvol*3*NTRIAL);
   else if (use_gpu==0)
      rnd=randgpu(nvol*3*NTRIAL);
   else
      rnd=NULL;

   return rnd;
}

double thermalize(double delta, int nsweep, int use_gpu)
{
   int i;
   double acc, *rnd;

   for (i=0; i<nsweep; i++)
   {
      rnd=rand4sweep(use_gpu);
      acc=metro_sweep(delta,rnd,use_gpu);
      delta=tune_delta(acc,delta);
   }

   return delta;
}

int main(int argc, char **argv)
{
   printf("%s Starting...\n", argv[0]);

   int i, nsweep,use_gpu;
   double acc, delta, s, iStart, mm, phi0, reh, s2;
   double *rnd;
   spin m;

   if (argc<MIN_NARG+1)
      usage();

   // read parameters from command line
   lambda=atof(argv[1]);
   kappa=atof(argv[2]);
   reh=atof(argv[3]);
   phi0=atof(argv[4]);
   nsweep=atoi(argv[5]);
   use_gpu=atoi(argv[6]);

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
   printf("gpu    = %d\n",use_gpu);
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
   set_eo();

   // allocate spins, set random values
   printf("Initalize spins... ");
   iStart=seconds();
   phi=(spin*)malloc(nvol*sizeof(spin));
   init_phi(phi0);
   printf("%f sec.\n\n",seconds()-iStart);

   if (use_gpu)
   {
      init_gpu();
   }

   // thermalize
   printf("Thermalize spins... ");
   iStart=seconds();
   delta=thermalize(delta,NTHERM,use_gpu);
   printf("%f sec.\n\n",seconds()-iStart);

   // Measurements
   printf("Measurements...\n\n");
   // Initialisierung des Statistikmoduls, vier Observablen
   clear5(4,500);

   s=action();
   m=magnet();
   mm=cuCabs(m)*cuCabs(m);

   printf("UPD\t A       \t S       \t RE(M)   \t IM(M)   \t |M|^2   \n");

   acc=0.0;
   s2=0.0;
   iStart=seconds();
   for (i=1; i<=nsweep; i++)
   {
      // Zufallszahlen
      rnd=rand4sweep(use_gpu);
      // sweep
      acc=metro_sweep(delta,rnd,use_gpu);
      // Messung
      s=action();
      m=magnet();
      mm=cuCabs(m)*cuCabs(m);
      printf("%d\t %f\t %f\t %f\t %f\t %f\n",i,acc,s,cuCreal(m),cuCimag(m),mm);
      // Messwerte für Auswertung speichern
      accum5(1,s);
      accum5(2,cuCreal(m));
      accum5(3,cuCimag(m));
      accum5(4,cuCabs(m)*cuCabs(m));

      s2+=s*s;
   }
   s2/=(double)nsweep;
   s2=sqrt((s2-aver5(1)*aver5(1))/(double)nsweep);

   printf("\n\n");
   printf("%d updates took %f sec. [sqrt(var(S)/N): %f, check: %f]\n\n",
                  nsweep,seconds()-iStart,s2,s2*sqrt(2.0*tauint5(1)));

   // Mittelwert, Fehler, Auto-Korrelationszeit
   printf("           Mittelwert, Fehler, integr. Auto-Korrelationszeit\n");
   printf(" S:        %f, %f, %f\n",aver5(1),sigma5(1),tauint5(1));
   printf(" RE(M):    %f, %f, %f\n",aver5(2),sigma5(2),tauint5(2));
   printf(" IM(M):    %f, %f, %f\n",aver5(3),sigma5(3),tauint5(3));
   printf(" |M|^2:    %f, %f, %f\n",aver5(4),sigma5(4),tauint5(4));

   free(lsize);
   free(nn[0]);
   free(nn);
   free(phi);
   randgpu(0);
}
