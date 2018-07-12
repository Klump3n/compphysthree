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
#include "added_stuff.h"
#include "added_stuff_gpu.h"

#include <cuda_runtime.h>

#define MIN_NARG 6

#define CHECK(call)                                             \
  {                                                             \
    const cudaError_t error = call;                             \
    if (error != cudaSuccess)                                   \
      {                                                         \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);  \
        fprintf(stderr, "code: %d, reason: %s\n", error,        \
                cudaGetErrorString(error));                     \
        exit(1);                                                \
      }                                                         \
  }

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

void gpu_stuff(int nsweep)
{

  // set up device
  int dev = 0;
  cudaDeviceProp deviceProp;
  CHECK(cudaGetDeviceProperties(&deviceProp, dev));
  printf("Using Device %d: %s\n", dev, deviceProp.name);
  CHECK(cudaSetDevice(dev));

  int *d_nn;
  CHECK(cudaMalloc((void**)&d_nn,nvol*(2*ndim+1)*sizeof(int)));
  CHECK(cudaMemcpy(d_nn, nn[0], nvol*(2*ndim+1)*sizeof(int), cudaMemcpyHostToDevice));

  CHECK(cudaMemcpyToSymbol(devLambda, &lambda, sizeof(double)));
  CHECK(cudaMemcpyToSymbol(devKappa, &kappa, sizeof(double)));
  CHECK(cudaMemcpyToSymbol(devNdim, &ndim, sizeof(int)));
  CHECK(cudaMemcpyToSymbol(devNvol, &nvol, sizeof(int)));

  int halfArrayLength = ceil((float) (nvol) / 2.);

  int *evenArray = (int *) calloc(halfArrayLength, sizeof(int));
  int *oddArray = (int *) calloc(halfArrayLength, sizeof(int));
  evenOddIndices(evenArray, oddArray);

  int *d_evenArray;
  CHECK(cudaMalloc((void**)&d_evenArray, halfArrayLength * sizeof(int)));
  CHECK(cudaMemcpy(d_evenArray, evenArray, halfArrayLength * sizeof(int), cudaMemcpyHostToDevice));
  int *d_oddArray;
  CHECK(cudaMalloc((void**)&d_oddArray, halfArrayLength * sizeof(int)));
  CHECK(cudaMemcpy(d_oddArray, oddArray, halfArrayLength * sizeof(int), cudaMemcpyHostToDevice));

  spin *d_bEvenArray;
  CHECK(cudaMalloc((void**)&d_bEvenArray, halfArrayLength * sizeof(spin)));
  CHECK(cudaMemset(d_bEvenArray, 0, halfArrayLength * sizeof(spin)));
  spin *d_bOddArray;
  CHECK(cudaMalloc((void**)&d_bOddArray, halfArrayLength * sizeof(spin)));
  CHECK(cudaMemset(d_bOddArray, 0, halfArrayLength * sizeof(spin)));

  spin *d_phi;
  CHECK(cudaMalloc((void**)&d_phi, nvol*sizeof(spin)));
  CHECK(cudaMemcpy(d_phi, phi, nvol*sizeof(spin), cudaMemcpyHostToDevice));

  int *d_accept;
  CHECK(cudaMalloc((void**)&d_accept, halfArrayLength * sizeof(int)));
  CHECK(cudaMemset(d_accept, 0, halfArrayLength * sizeof(int)));

  spin *d_phi_intermediate;
  CHECK(cudaMalloc((void**)&d_phi_intermediate, nvol*sizeof(spin)));
  CHECK(cudaMemset(d_phi_intermediate, 0, nvol*sizeof(spin)));

  double *d_aloc_comp;
  CHECK(cudaMalloc((void**)&d_aloc_comp, halfArrayLength * sizeof(double)));
  CHECK(cudaMemset(d_aloc_comp, 0, halfArrayLength * sizeof(double)));

  double *d_aloc_calc;
  CHECK(cudaMalloc((void**)&d_aloc_calc, halfArrayLength * sizeof(double)));
  CHECK(cudaMemset(d_aloc_calc, 0, halfArrayLength * sizeof(double)));

  /* NTRIAL IS 10!!! */
  double *d_rnd = randgpu_device_ptr(3*nvol*10);

  double *rnd = (double *) calloc(nvol, sizeof(double));
  CHECK(cudaMemcpy(rnd, d_rnd, nvol*3*10*sizeof(double), cudaMemcpyDeviceToHost));
  /* for (int i=0; i<nvol*30; i++) */
  /*   { */
  /*     printf("%d, %f\n", i, doublearray[i]); */
  /*   } */

  spin m, gpu_m;
  double mm, gpu_mm;

  double delta = .2;
  double cpu_delta = .2;

  double acc = 0.0;
  double cpu_acc = 0.0;

  spin *backup_phi = (spin *) malloc(nvol*sizeof(spin));

  printf("ITER \t ACC \t\t ACC_C \t\t DELTA \t\t DELTA_C \t GPU_MM \t CPU_MM\n");
  for (int i=1; i<=nsweep; i++)
    {

//      cpu_acc = metro_sweep_alt(cpu_delta, evenArray, oddArray, rnd);
//      cpu_delta=tune_delta(cpu_acc,cpu_delta);
//
//      m=magnet();
//      mm=cuCabs(m)*cuCabs(m);

      acc=gpu_sweep(d_phi,
                    d_evenArray,
                    d_oddArray,
                    d_bEvenArray,
                    d_bOddArray,
                    d_nn,
                    d_accept,
                    d_phi_intermediate,
                    d_aloc_comp,
                    d_aloc_calc,
                    d_rnd,
                    delta
                    );

//      memcpy(backup_phi, phi, nvol*sizeof(spin)); /* keep original phi */
//      CHECK(cudaMemcpy(phi, d_phi, nvol*sizeof(spin), cudaMemcpyDeviceToHost));
      gpu_m = magnet();
//      memcpy(phi, backup_phi, nvol*sizeof(spin)); /* restore phi */

      gpu_mm=cuCabs(gpu_m)*cuCabs(gpu_m);

      delta=tune_delta(acc,delta);

      /* printf("%d\t %f\t %f\t %f\t %f\t %f\t %f\n",i,acc,delta,s,cuCreal(m),cuCimag(m),mm); */

      printf("%d \t %f \t %f \t %f \t %f \t %f \t %f \n", i, acc, cpu_acc, delta, cpu_delta, gpu_mm, mm);
    }

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
   CHECK(cudaMemcpyToSymbol(devH, &h, sizeof(spin)));

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


   spin *initial_phi = (spin *) malloc(nvol*sizeof(spin));
   memcpy(initial_phi, phi, nvol*sizeof(spin)); /* keep original phi */
   /* printf("i phi init_phi\n"); */
   /* for (int i=0; i<nvol; i++) */
   /*   { */
   /*     printf("%d, %f, %f\n", i, cuCabs(phi[i]), cuCabs(initial_phi[i])); */
   /*   } */


   /* s=action(); */
   /* m=magnet(); */
   /* mm=cuCabs(m)*cuCabs(m); */

   /* printf("UPD\t A       \t DELTA   \t S       \t RE(M)   \t IM(M)   \t |M|^2   \n"); */
   /* printf("%d\t %f\t %f\t %f\t %f\t %f\t %f\n",0,0.0,delta,s,cuCreal(m),cuCimag(m),mm); */

   /* acc=0.0; */
   /* iStart=seconds(); */
   /* for (i=1; i<=nsweep; i++) */
   /* { */
   /*    acc=metro_sweep(delta); */
   /*    delta=tune_delta(acc,delta); */
   /*    s=action(); */
   /*    m=magnet(); */
   /*    mm=cuCabs(m)*cuCabs(m); */
   /*    printf("%d\t %f\t %f\t %f\t %f\t %f\t %f\n",i,acc,delta,s,cuCreal(m),cuCimag(m),mm); */
   /* } */

   /* printf("\n\n"); */
   /* printf("%d updates took %f sec.\n\n",nsweep,seconds()-iStart); */



   memcpy(phi, initial_phi, nvol*sizeof(spin)); /* restore phi */
   /* printf("i phi init_phi\n"); */
   /* for (int i=0; i<nvol; i++) */
   /*   { */
   /*     printf("%d, %f, %f\n", i, cuCabs(phi[i]), cuCabs(initial_phi[i])); */
   /*   } */
   gpu_stuff(nsweep);





   free(lsize);
   free(nn[0]);
   free(nn);
   free(phi);
}
