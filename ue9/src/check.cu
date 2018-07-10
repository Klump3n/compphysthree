#define DEFINE_GLOBAL

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "global.h"
#include "geom_pbc.h"
#include "spin.h"
#include "randgpu.h"
#include "metropolis.h"

#include "added_stuff.h"

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

int check_action1(double a, double b)
{
   int idx;
   double act1,act2,tmp;

   for (idx=0; idx<nvol; idx++)
   {
      phi[idx]=make_spin(a,b);
      //printf("%f + %f * i\n", cuCreal(phi[idx]), cuCimag(phi[idx]));
   }

   act1=action();
   tmp=a*a+b*b;
   act2=((double)nvol)*(   (1.0-2*lambda-2.0*kappa*ndim)*tmp
                         + lambda*(1+tmp*tmp)
                         - 2.0*(a*cuCreal(h)+b*cuCimag(h))  );


   printf(" Diff: %e (%e)\n",fabs((act1-act2)/act2),(sqrt(nvol*(100+50*ndim))*DBL_EPSILON));

   return (fabs((act1-act2)/act2)<(nvol*DBL_EPSILON));
}

int check_action2(double alpha)
{
   int idx;
   double act1,act2;
   spin f;

   random_cnfg();

   f=make_spin(cos(alpha),sin(alpha));
   h=make_spin(0.0,0.0);

   act1=action();

   for (idx=0; idx<nvol; idx++)
   {
      phi[idx]=cuCmul(phi[idx],f);
   }

   act2=action();

   printf(" Diff: %e\n",fabs((act1-act2)/act2));

   return (fabs((act1-act2)/act2)<sqrt(nvol)*DBL_EPSILON);
}

int check_magnet(double alpha)
{
   int idx;
   double d;
   spin m1,m2,f;

   random_cnfg();

   f=make_spin(cos(alpha),sin(alpha));

   m1=magnet();
   m1=cuCmul(m1,f);

   for (idx=0; idx<nvol; idx++)
   {
      phi[idx]=cuCmul(phi[idx],f);
   }
   m2=magnet();

   d=cuCabs(cuCsub(m1,m2))/cuCabs(m2);

   printf(" Diff: %e\n",d);

   return (d<sqrt(nvol)*DBL_EPSILON);
}

int check_alocal()
{
   int idx,ifail;
   double act1,act2,a1,a2,a3,diff,mdiff;
   spin tmp;

   random_cnfg();

   ifail=0;
   mdiff=0.0;
   act1=action();
   for (idx=0; idx<nvol; idx++)
   {
      a1=alocal(idx);
      a3=alocal2(idx);
      if (fabs(a1-a3)!=0.0)
      {
         printf("idx: %d, diff(a1-a3): %e\n",idx,fabs(a1-a3));
         ifail=1;
      }

      tmp=phi[idx];
      phi[idx]=make_spin(cuCreal(phi[idx])+10.0,cuCimag(phi[idx]));
      act2=action();
      a2=alocal(idx);
      phi[idx]=tmp;

      diff=fabs(((-act2+act1)-(-a2+a1))/(-act2+act1));
      if (diff>1e-7)
      {
         printf("idx: %d, diff: %e\n",idx,diff);
         ifail=1;
      }
      if (diff>mdiff)
         mdiff=diff;
   }

   printf(" Max. diff: %e\n",mdiff);

   return (ifail==0);
}

void added_stuff_checks()
{
  int *target_array = (int *) malloc(ndim*sizeof(int));

  int index = 0;
<<<<<<< HEAD
  int evenOdd = 0;
=======
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
  memset(target_array, 0, ndim*sizeof(int));
  indexArray(index, target_array);
  printf("%d: ", index);
  for (int j=1; j<=ndim; j++)
    {
      printf("%d ", target_array[j]);
    }
<<<<<<< HEAD
  evenOdd = indexEvenOdd(target_array);
  if (evenOdd == 0)
    {
      printf("even");
    }
  else
    {
      printf("odd");
    }
=======
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
  printf("\n");

  index = 1;
  memset(target_array, 0, ndim*sizeof(int));
  indexArray(index, target_array);
  printf("%d: ", index);
    for (int j=1; j<=ndim; j++)
      {
        printf("%d ", target_array[j]);
      }
<<<<<<< HEAD
    evenOdd = indexEvenOdd(target_array);
    if (evenOdd == 0)
      {
        printf("even");
      }
    else
      {
        printf("odd");
      }
=======
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
  printf("\n");

  index = 7;
  memset(target_array, 0, ndim*sizeof(int));
  indexArray(index, target_array);
  printf("%d: ", index);
    for (int j=1; j<=ndim; j++)
      {
        printf("%d ", target_array[j]);
      }
<<<<<<< HEAD
    evenOdd = indexEvenOdd(target_array);
    if (evenOdd == 0)
      {
        printf("even");
      }
    else
      {
        printf("odd");
      }
=======
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
  printf("\n");

  index = 8;
  memset(target_array, 0, ndim*sizeof(int));
  indexArray(index, target_array);
  printf("%d: ", index);
    for (int j=1; j<=ndim; j++)
      {
        printf("%d ", target_array[j]);
      }
<<<<<<< HEAD
    evenOdd = indexEvenOdd(target_array);
    if (evenOdd == 0)
      {
        printf("even");
      }
    else
      {
        printf("odd");
      }
=======
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
  printf("\n");

  index = 9;
  memset(target_array, 0, ndim*sizeof(int));
  indexArray(index, target_array);
  printf("%d: ", index);
    for (int j=1; j<=ndim; j++)
      {
        printf("%d ", target_array[j]);
      }
<<<<<<< HEAD
    evenOdd = indexEvenOdd(target_array);
    if (evenOdd == 0)
      {
        printf("even");
      }
    else
      {
        printf("odd");
      }
=======
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
  printf("\n");

  index = 24;
  memset(target_array, 0, ndim*sizeof(int));
  indexArray(index, target_array);
  printf("%d: ", index);
    for (int j=1; j<=ndim; j++)
      {
        printf("%d ", target_array[j]);
      }
<<<<<<< HEAD
    evenOdd = indexEvenOdd(target_array);
    if (evenOdd == 0)
      {
        printf("even");
      }
    else
      {
        printf("odd");
      }
=======
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
  printf("\n");

  index = 63;
  memset(target_array, 0, ndim*sizeof(int));
  indexArray(index, target_array);
  printf("%d: ", index);
  for (int j=1; j<=ndim; j++)
    {
      printf("%d ", target_array[j]);
    }
<<<<<<< HEAD
  evenOdd = indexEvenOdd(target_array);
  if (evenOdd == 0)
    {
      printf("even");
    }
  else
    {
      printf("odd");
    }
=======
>>>>>>> ddc12ebc74ec86c0fa872a0f6877224586b41ed4
  printf("\n");

}

int check_update_spin(double delta, double rnd2)
{
   int idx,acc,ifail;
   double *rnd,diff,mdiff=0.0;
   spin tmp,tmp2;

   ifail=0;
   for (idx=0; idx<nvol; idx++)
   {
   	rnd=randgpu(3);
   	rnd[2]=rnd2; // =0.0: accept in any case
   	acc=0;
   	tmp=phi[idx];
   	acc+=update_spin(idx,rnd,delta,1);
      if (acc==1)
   	{
         tmp2=cuCsub(tmp,phi[idx]);
         tmp2=cuCadd(tmp2,make_spin(delta*(rnd[0]-0.5),delta*(rnd[1]-0.5)));
         diff=cuCabs(tmp2);
         if (diff>sqrt(6)*DBL_EPSILON)
         {
   	      printf("idx: %d, diff1: %e\n",idx,diff);
   	      ifail=1;
            if (diff>mdiff)
               mdiff=diff;
         }
         // flip delta
         rnd[2]=0.0; // =0.0: accept in any case next update
      	acc+=update_spin(idx,rnd,-delta,1);
   	}
   	if (acc!=2 && !(rnd2!=0.0 && acc==0))
   	{
   	   printf("idx: %d, acc: %d\n",idx,acc);
   	   ifail=1;
   	}
   	diff=cuCabs(cuCsub(tmp,phi[idx]))/cuCabs(tmp);
      if (diff>sqrt(6)*DBL_EPSILON)
	   {
         printf("idx: %d, diff2: %e\n",idx,diff);
	      ifail=1;
	   }
      if (diff>mdiff)
         mdiff=diff;
   }

   printf(" Max. diff: %e\n",mdiff);

   return (ifail==0);
}

int main(int argc, char **argv)
{
   printf("%s Starting...\n", argv[0]);

   int i;
   double alpha;

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

   if (ndim==2 && nvol<256)
   {
      print_nn(0);
      print_nn(1);
      print_nn(2);
      print_nn(3);
      print_nn(4);
   }

   phi=(spin*)malloc(nvol*sizeof(spin));

   printf("Action:\n");
   h=make_spin(0.3,0.5);
   lambda=0.7;
   kappa=0.06;
   if (check_action1(0.5,0.0))
      printf(" Check1 erfolgreich.\n");
   else
      printf(" Check1 fehlgeschalgen!.\n");
   alpha=0.45;
   if (check_action2(alpha))
      printf(" Check2 erfolgreich.\n\n");
   else
      printf(" Check2 fehlgeschalgen!.\n\n");

   printf("Magnetisierung:\n");
   alpha=0.45;
   if (check_magnet(alpha))
      printf(" Check erfolgreich.\n\n");
   else
      printf(" Check fehlgeschalgen!.\n\n");

   printf("loakle Verteilung/Wirkung:\n");
   if (check_alocal())
      printf(" Check erfolgreich.\n\n");
   else
      printf(" Check fehlgeschalgen!.\n\n");

   printf("Updates des Spins:\n");
   if (check_update_spin(0.0,0.0)) // zero delta, accept always
      printf(" Check1 erfolgreich.\n");
   else
      printf(" Check1 fehlgeschalgen!.\n");
   if (check_update_spin(0.1,0.0)) // non-zero delta, accept always
      printf(" Check2 erfolgreich.\n");
   else
      printf(" Check2 fehlgeschalgen!.\n");
   if (check_update_spin(0.1,0.5)) // non-zero delta, accept sometimes
      printf(" Check3 erfolgreich.\n");
   else
      printf(" Check3 fehlgeschalgen!.\n");
   if (check_update_spin(0.1,1000000.0)) // non-zero delta, never accept
      printf(" Check4 erfolgreich.\n\n");
   else
      printf(" Check4 fehlgeschalgen!.\n\n");

   printf("Checking added stuff\n");
   added_stuff_checks();

   free(lsize);
   free(nn[0]);
   free(nn);
   free(phi);
}
