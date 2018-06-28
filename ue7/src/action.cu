#define DEFINE_GLOBAL

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "global.h"
#include "geom_pbc.h"
#include "randgpu.h"

#include "tests.h"

void print_nn(int k)
{
   int i,j,idx;

   if (k>0 && k<3)
      printf("Naechste Nachbarn, %d-Richtung:\n",k);
   else if (k>2 && k<5)
      printf("Naechste Nachbarn, -%d-Richtung:\n",k-2);
   else
      printf("Fortlaufender Index:\n");
   idx=0;
   for (j=0; j<lsize[2]; j++)
   {
      printf("  ");
      for (i=0; i<lsize[1]; i++)
      {
         if (k>0)
            printf("%d ",nn[k][idx]);
         else
            printf("%d ",idx);
         idx+=1;
      }
      printf("\n");
   }
}

double action(double lambda, double kappa, cuDoubleComplex h)
{
   int idx,k;
   double act,tmp;

   act=0.0;
   kappa*=2.0;
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
      act-=kappa*tmp;
      act-=2.0*(cuCreal(phi[idx])*cuCreal(h) + cuCimag(phi[idx])*cuCimag(h));
   }

   return act;
}

double alocal(int idx, double lambda, double kappa, cuDoubleComplex h)
{
   int k;
   double a,tmp;
   cuDoubleComplex b,tmpc;

   b=h;
   for (k=1; k<=ndim; k++)
   {
      tmpc=cuCadd(phi[nn[k][idx]],phi[nn[ndim+k][idx]]);
      b=make_cuDoubleComplex(cuCreal(b)+kappa*cuCreal(tmpc),cuCimag(b)+kappa*cuCimag(tmpc));
   }

   tmp=cuCreal(phi[idx])*cuCreal(phi[idx]) + cuCimag(phi[idx])*cuCimag(phi[idx]);
   a=2.0*(cuCreal(b)*cuCreal(phi[idx])+cuCimag(b)*cuCimag(phi[idx]))-tmp;
   tmp-=1.0;
   a-=lambda*tmp*tmp;

   #ifdef DEBUG
      printf("b: %f, a: %f\n",b,a);
   #endif

   return -a;
}

int check1(double lambda, double kappa, cuDoubleComplex h, double a, double b)
{
   int idx;
   double act1,act2,tmp;

   for (idx=0; idx<nvol; idx++)
   {
      phi[idx]=make_cuDoubleComplex(a,b);
      //printf("%f + %f * i\n", cuCreal(phi[idx]), cuCimag(phi[idx]));
   }

   act1=action(lambda,kappa,h);
   tmp=a*a+b*b;
   act2=((double)nvol)*(   (1.0-2*lambda-2.0*kappa*ndim)*tmp
                         + lambda*(1+tmp*tmp)
                         - 2.0*(a*cuCreal(h)+b*cuCimag(h))  );


   printf("Check1: %e (%e)\n",fabs((act1-act2)/act2),(sqrt(nvol*(100+50*ndim))*DBL_EPSILON));

   return (fabs((act1-act2)/act2)<(nvol*DBL_EPSILON));
}

void random_cnfg(void)
{
   int idx;

   for (idx=0; idx<nvol; idx++)
   {
      phi[idx]=make_cuDoubleComplex((double)(rand() & 0xFF ) / 99.0,(double)(rand() & 0xFF ) / 99.0);
   }
}

int check2(double lambda, double kappa, double alpha)
{
   int idx;
   double act1,act2;
   cuDoubleComplex f,h;

   random_cnfg();

   f=make_cuDoubleComplex(cos(alpha),sin(alpha));
   h=make_cuDoubleComplex(0.0,0.0);

   act1=action(lambda,kappa,h);

   for (idx=0; idx<nvol; idx++)
   {
      phi[idx]=cuCmul(phi[idx],f);
   }

   act2=action(lambda,kappa,h);

   printf("Check2: %e\n",fabs((act1-act2)/act2));

   return (fabs((act1-act2)/act2)<sqrt(nvol)*DBL_EPSILON);
}

int check_alocal(double lambda, double kappa, cuDoubleComplex h)
{
   int idx,ifail;
   double act1,act2,a1,a2,diff,mdiff;
   cuDoubleComplex tmp;

   random_cnfg();

   ifail=0;
   mdiff=0.0;
   act1=action(lambda,kappa,h);
   for (idx=0; idx<nvol; idx++)
   {
      tmp=phi[idx];
      a1=alocal(idx,lambda,kappa,h);
      phi[idx]=make_cuDoubleComplex(cuCreal(phi[idx])+10.0,cuCimag(phi[idx]));
      act2=action(lambda,kappa,h);
      a2=alocal(idx,lambda,kappa,h);
      phi[idx]=tmp;

      diff=fabs(((-act2+act1)-(-a2+a1))/(-act2+act1));
      if (diff>1e-7)
      {
         printf("idx: %d, diff: %e %e\n",idx,diff,(-act2+act1)-(-a2+a1));
         ifail=1;
      }
      if (diff>mdiff)
         mdiff=diff;
   }

   printf("Check alocal: Max. diff: %e\n",mdiff);

   return (ifail==0);
}

int main(int argc, char **argv)
{
   printf("%s Starting...\n", argv[0]);

   int i;
   double lambda, kappa, alpha;
   cuDoubleComplex h;

   if (argc>1)
   {
      ndim=argc-1;
      lsize=(int*)malloc((ndim+1)*sizeof(int));
      for (i=1; i<argc; i++)
      {
         lsize[i]=atoi(argv[i]);
      }
   }
   else
   {
      ndim=1;
      lsize=(int*)malloc((ndim+1)*sizeof(int));
      lsize[1]=8;
   }

   printf("Gittergroesse: %d",lsize[1]);
   for (i=2; i<=ndim; i++)
   {
      printf(" x %d",lsize[i]);
   }
   printf("\n\n");

   geom_pbc();

   if (ndim==2)
   {
      print_nn(0);
      print_nn(1);
      print_nn(2);
      print_nn(3);
      print_nn(4);
   }

   phi=(cuDoubleComplex*)malloc(nvol*sizeof(cuDoubleComplex));

   h=make_cuDoubleComplex(0.3,0.5);
   lambda=0.7;
   kappa=0.06;
   if (check1(lambda,kappa,h,0.5,0.0))
      printf("Check1 erfolgreich.\n");
   else
      printf("Check1 fehlgeschalgen!.\n");

   alpha=0.45;
   if (check2(lambda,kappa,alpha))
      printf("Check2 erfolgreich.\n");
   else
      printf("Check2 fehlgeschalgen!.\n");

   if (check_alocal(lambda,kappa,h))
      printf("Check alocal erfolgreich.\n");
   else
      printf("Check alocal fehlgeschalgen!.\n");

   printf("\n");

   double *rnd;
   rnd=randgpu(20);
   for (int i=1; i<20; i++)
   {
      printf(" %.6f\n",rnd[i]);
   }

//   mag_test();
   /* other_test(); */
   /* int fitting_res = delta_fitting_test(); */
//   int set_res = spin_set_test();
//   spin_update_test();
  boltzmag();

   free(lsize);
   free(nn[0]);
   free(nn);
   free(phi);
}
