#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
   Globale Variablen stehen in allen Funktionen zur Verfuegung.
   Achtung: Das gilt *nicht* fuer Kernel-Funktionen!
*/
int Nx, Ny, npts;
int *active;

/*
   Fuer die Koordinaten:
      i = 0,1,...,Nx+1
      j = 0,1,...,Ny+1
   wird der fortlaufenden Index berechnet
*/
int coord2index(int i, int j)
{
   return j*(Nx+2) + i;
}

/*
   Das Flag-Array der aktiven/inneren Punkte wird gesetzt.
*/
void active_pts()
{
   int idx,i,j;

   active=(int*)malloc(npts*sizeof(int));

   idx=0; // fortlaufender Index
   for (j=0; j<Ny+2; j++)
   {
      for (i=0; i<Nx+2; i++)
      {
         if ((i==0)||(j==0)||(i==Nx+1)||(j==Ny+1))
            active[idx]=0; // Randpunkt
         else
            active[idx]=1; // innerer Punkt

         idx+=1;
      }
   }
}

/*
   Der Vektor p wird im Inneren auf zufaellige Werte gesetzt
*/
void random_vector(double *p)
{
   int idx;

   for(idx = 0; idx < npts; idx++)
   {
      if (active[idx])
         p[idx] = (double)(rand() & 0xFF ) / 10.0;
   }
}

/*
   Das Flag-Array der aktiven/inneren Punkte wird als
   2D Gitter ausgegeben.
*/
void print_active()
{
   int i,j,idx;

   printf("active points:\n");
   idx=0;
   for (j=0; j<Ny+2; j++)
   {
      printf("  ");
      for (i=0; i<Nx+2; i++)
      {
         printf("%d ",active[idx]);
         idx+=1;
      }
      printf("\n");
   }
}


/*
   Norm-Quadrat vom Vektor v.
*/
double norm_sqr(double *v)
{
   int idx;
   double r=0.0;
   for (idx=0; idx<npts; idx++)
   {
      r+=v[idx]*v[idx];
   }
   return r;
}

/*
   Der Vektor p wird als 2D Gitter fuer i,j<=16 ausgegeben. Es werden innere/aktive
   und, falls flag>0, auch die aeusseren Punkte ausgegeben.
*/
void print_vector(char *name, double *p, int flag)
{
   int i,j,idx;
   double nrm;

   printf("%s = \n",name);
   idx=0;
   for (j=0; j<Ny+2; j++)
   {
      if (j>16)
      {
         printf("  ...\n");
         break;
      }
      printf("  ");
      for (i=0; i<Nx+2; i++)
      {
         if ((i<16)&&((flag>0)||(active[idx])))
           printf("%.2f ",p[idx]);
         if (i==16)
           printf("...");
         idx+=1;
      }
      printf("\n");
   }
   nrm=norm_sqr(p);
   printf("||%s|| = %.8f\n",name,sqrt(nrm));
}

/*
  Die geforderte Funktion
 */
void laplace_2d(double *w, double *v) {

  int i, j, n = Nx+2;
  int npts = n * n;

  /* two loops over every element of the unknowns */
  for (i=0; i<npts; i++) {
    for (j=0; j<npts; j++) {

      /* skip the boundary */
      if (
          (i % n == 0) ||
          ((int) (i / n) == 0) ||
          (i % n == n-1) ||
          ((int) (i / n) == n-1)
          ) {
        continue;
      }

      /* diagonal */
      if (i == j) {
        w[i] += -4 * v[j];
      }

      /* first two off diagonals */
      if (
          ((i == j+1) && (i % n != 0))
          ||
          ((i == j-1) && (j % n != 0))
          )
        {
          w[i] += 1 * v[j];
        }

      /* other two off diagonals */
      if ((i == j+n) || (i == j-n)) {
        w[i] += 1 * v[j];
      }
    }
  }
}

double scalar_product(double *u, double *v) {
  int n = (Nx + 2)*(Ny + 2);
  int i;

  double w = 0;

  for (i=0; i<npts; i++) {
    printf("%d, %f, %f, %f\n", i, u, v, w);
    w += u[i] * v[i];
  }

  return w;
}

void vector_addition(double *v,double *u, double *w) {
    int i;
    for (i=0; i<npts; i++) {
        w[i] += u[i] + v[i];
    }      
}

void scale_vector(double prefactor, double *v, double *w) {
  int n = (Nx + 2)*(Ny + 2);
  int i;

  for (i=0; i<npts; i++) {
    w[i] = prefactor * v[i];
  }
}

void conjugate_gradient_laplace(double *x, double *b) {
  /* int i = 0; */

  double tol = 1e-15;
  double alpha, beta;

  int nBytes=npts*sizeof(double);

  /* double* x =(double*)malloc(npts*sizeof(double)); /\* residuum *\/ */
  double* r =(double*)malloc(npts*sizeof(double)); /* residuum */
  double* p =(double*)malloc(npts*sizeof(double)); /* residuum */
  double* Ap =(double*)malloc(npts*sizeof(double)); /* residuum */

  double rr_old, rr_new;

  /* memset(x, 0, nBytes); */
  memset(r, 0, nBytes);
  memset(p, 0, nBytes);
  memset(Ap, 0, nBytes);

  r = b;
  print_vector("r",r,1);
  /* *p = *r; */

  printf("%f\n", norm_sqr(r));
  rr_old = scalar_product(r, r);

  int j;
  /* for (j=0; *rr_old<tol; j++) { */
  for (j=0; j<1000; j++) {

    /* /\* memset(x, 0, nBytes); *\/ */
      /* memset(r, 0, nBytes); */
      /* memset(p, 0, nBytes); */
      /* memset(Ap, 0, nBytes); */
      /* printf("%d\n", j); */
      /* abbruch bed */

    if (j > 100) {
        printf("iter > 100\n");

        break;
      }
  }

  /* free(x); */
  free(r);
  free(p);
  free(Ap);
}

int main(int argc, char **argv)
{
   printf("%s Starting...\n", argv[0]);

   int nBytes;
   double *w, *v;

   // Globale Variablen setzen:
   // Anzahl der Inneren Punkte in x- und y-Richtung
   Nx=8;
   Ny=8;
   // Gesamtanzahl der Gitterpunkte
   npts=(Nx+2)*(Ny+2);
   // Aktive Punkte - Array
   active_pts();

   // Speicherbedarf pro Vektor in Byte
   nBytes=npts*sizeof(double);

   // Speicher für Vektoren allozieren
   w=(double*)malloc(npts*sizeof(double));
   v=(double*)malloc(npts*sizeof(double));

   double* x;
   x=(double*)malloc(npts*sizeof(double));
   memset(x, 0, nBytes);

   // auf Null setzen
   memset(w, 0, nBytes);
   memset(v, 0, nBytes);

   // Aktive Punkte ausgeben
   if ((Nx<=16)&&(Ny<=16))
      print_active();

   // Einheitsvektor
   v[coord2index(Nx/2+1,Nx/2+1)]=1.0; // v=0, ausser am Gitterpunkt (Nx/2+1,Ny/2+1)
   print_vector("v",v,1);

   // Zufaelliger Vektor
   random_vector(v);
   print_vector("v",v,1);

   laplace_2d(w, v);
   print_vector("w", w, 1);

   /* conjugate_gradient_laplace(x, v); */

   free(active);
   free(w);
   free(v);

   return (0);
}
