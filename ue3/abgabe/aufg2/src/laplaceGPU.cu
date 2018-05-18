#include <cuda_runtime.h>
#include <stdio.h>
#include <sys/time.h>


/*
  Globale Variablen stehen in allen Funktionen zur Verfuegung.
  Achtung: Das gilt *nicht* fuer Kernel-Funktionen!
*/
int Nx, Ny, N, npts;
int *active;

/*
 * Dieses Beispiel demonstriert die Addition zweier Arrays.
 * addArrayGPU soll die Arbeit ueber CUDA Threads auf der GPU verteilen.
 * addArrayHost iteriert sequentiell durch die Vektorelemente auf dem Host.
 */

// Macro zur Fehlerauswertung
#define CHECK(call)                                                            \
{                                                                              \
    const cudaError_t error = call;                                            \
    if (error != cudaSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                cudaGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}

double seconds()
{
  struct timeval tp;
  struct timezone tzp;
  int i = gettimeofday(&tp, &tzp);
  return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

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
void random_vector(float *p)
{
   int idx;

   for(idx = 0; idx < npts; idx++)
   {
      if (active[idx])
         p[idx] = (float)(rand() & 0xFF ) / 10.0;
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
float norm_sqr(float *v)
{
   int idx;
   float r=0.0;
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
void print_vector(char *name, float *p, int flag)
{
   int i,j,idx;
   float nrm;

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

__global__
void laplace2d_GPU(float *w, float *v, const int N) {

  /* calculate the id */
  int blockOffset = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y;
  int threadId = threadIdx.x + threadIdx.y * blockDim.x;
  int idx = blockOffset + threadId;

  int first_coord = (int) (idx % N);
  int second_coord = (int) (idx / N);

  if (
      (first_coord > 0) && (first_coord < (N-1)) &&
      (second_coord > 0) && (second_coord < (N-1))
      ) {
    w[idx] = -4. * v[idx] + v[idx-1] + v[idx+1] + v[idx-N] + v[idx+N];
    /* printf("%d, %f, %f\n", idx, w[idx], v[idx]); */
  }
}

/*
 * laplace_2d(float *w, float *v)
 *
 * Calculates the product A*v for a vector v and writes the result into w.
 * A is a laplacian.
 *
 */
void laplace2d_CPU(float *w, float *v) {

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

/*
 * vector_addition(float *v, float *u, float *w)
 *
 * Adds two vectors u and v together and writes the result into w.
 *
 */
void vector_addition_CPU(float *u, float *v, float *w) {
  int i;

  for (i=0; i<npts; i++) {
    w[i] += u[i] + v[i];
  }
}

__global__
void vector_addition_GPU(float *u, float *v, float *w) {

  /* calculate the id */
  int blockOffset = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y;
  int threadId = threadIdx.x + threadIdx.y * blockDim.x;
  int idx = blockOffset + threadId;

  w[idx] = u[idx] + v[idx];
}

/*
 * scale_vector(float prefactor, float *v, float *w)
 *
 * Scales a vector v by a prefactor and writes the result into w.
 *
 */
void scale_vector_CPU(float prefactor, float *v, float *w) {
  int i;

  for (i=0; i<npts; i++) {
    w[i] = prefactor * v[i];
  }
}

__global__
void scale_vector_GPU(float prefactor, float *v, float *w) {

  /* calculate the id */
  int blockOffset = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y;
  int threadId = threadIdx.x + threadIdx.y * blockDim.x;
  int idx = blockOffset + threadId;

  w[idx] = prefactor * v[idx];
}

void laplaceOneLoop(FILE* laplace_speedup_file, const int N, int gridX, int gridY, int threadX, int threadY) {

  int nBytes;
  float *h_w, *h_v, *h_u;

  // Globale Variablen setzen:
  // Anzahl der Inneren Punkte in x- und y-Richtung
  Nx = N;
  Ny = N;

  // Gesamtanzahl der Gitterpunkte
  npts=(Nx+2)*(Ny+2);
  // Aktive Punkte - Array
  active_pts();

  // Speicherbedarf pro Vektor in Byte
  nBytes=npts*sizeof(float);

  // Speicher für Vektoren allozieren
  h_w = (float *) malloc(npts * sizeof(float));
  h_v = (float *) malloc(npts * sizeof(float));
  h_u = (float *) malloc(npts * sizeof(float));

  // auf Null setzen
  memset(h_w, 0, nBytes);
  memset(h_v, 0, nBytes);
  memset(h_u, 0, nBytes);

  // Zufaelliger Vektor
  random_vector(h_v);

  /* print_vector("v",h_v,1); */

  // Device-Speicher allozieren mit cudaMalloc
  float *d_v, *d_w;
  CHECK(cudaMalloc((float**)&d_v, nBytes));
  CHECK(cudaMalloc((float**)&d_w, nBytes));

  // kopieren Host -> Device mit cudaMemcpy
  CHECK(cudaMemcpy(d_v, h_v, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_w, h_w, nBytes, cudaMemcpyHostToDevice));

  dim3 block(gridX, gridY);
  dim3 grid(threadX, threadY);

  double t_GPU_start = seconds();
  laplace2d_GPU <<<block,grid>>> (d_w, d_v, N+2);
  cudaDeviceSynchronize();
  double t_GPU_end = seconds();
  double t_GPU = t_GPU_end - t_GPU_start;

  /* kopieren Device -> Host mit cudaMemcpy */
  CHECK(cudaMemcpy(h_w, d_w, nBytes, cudaMemcpyDeviceToHost));

  /* print_vector("w_GPU",h_w,1); */

  // Device-Speicher freigeben
  CHECK(cudaFree(d_v));
  CHECK(cudaFree(d_w));

  double t_CPU_start = seconds();
  laplace2d_CPU(h_u, h_v);
  double t_CPU_end = seconds();
  double t_CPU = t_CPU_end - t_CPU_start;

  /* print_vector("w_CPU",h_w,1); */

  fprintf(
          laplace_speedup_file,
          "%lf, %lf, %lf,\n",
          t_GPU, t_CPU, t_CPU/t_GPU
          );

  free(active);
  free(h_w);
  free(h_v);

}

void vectorScaleOneLoop(FILE* vector_scale_speedup_file, const int N, int gridX, int gridY, int threadX, int threadY) {

  int nBytes;
  float *h_w, *h_v, *h_u;

  // Globale Variablen setzen:
  // Anzahl der Inneren Punkte in x- und y-Richtung
  Nx = N;
  Ny = N;

  // Gesamtanzahl der Gitterpunkte
  npts=(Nx+2)*(Ny+2);
  // Aktive Punkte - Array
  active_pts();

  // Speicherbedarf pro Vektor in Byte
  nBytes=npts*sizeof(float);

  // Speicher für Vektoren allozieren
  h_w = (float *) malloc(npts * sizeof(float));
  h_v = (float *) malloc(npts * sizeof(float));
  h_u = (float *) malloc(npts * sizeof(float));

  // auf Null setzen
  memset(h_w, 0, nBytes);
  memset(h_v, 0, nBytes);
  memset(h_u, 0, nBytes);

  // Zufaelliger Vektor
  random_vector(h_v);

  /* print_vector("v",h_v,1); */

  // Device-Speicher allozieren mit cudaMalloc
  float *d_v, *d_w;
  CHECK(cudaMalloc((float**)&d_v, nBytes));
  CHECK(cudaMalloc((float**)&d_w, nBytes));

  // kopieren Host -> Device mit cudaMemcpy
  CHECK(cudaMemcpy(d_v, h_v, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_w, h_w, nBytes, cudaMemcpyHostToDevice));

  dim3 block(gridX, gridY);
  dim3 grid(threadX, threadY);

  double t_GPU_start = seconds();
  laplace2d_GPU <<<block,grid>>> (d_w, d_v, N+2);
  cudaDeviceSynchronize();
  double t_GPU_end = seconds();
  double t_GPU = t_GPU_end - t_GPU_start;

  /* kopieren Device -> Host mit cudaMemcpy */
  CHECK(cudaMemcpy(h_w, d_w, nBytes, cudaMemcpyDeviceToHost));

  /* print_vector("w_GPU",h_w,1); */

  // Device-Speicher freigeben
  CHECK(cudaFree(d_v));
  CHECK(cudaFree(d_w));

  double t_CPU_start = seconds();
  laplace2d_CPU(h_u, h_v);
  double t_CPU_end = seconds();
  double t_CPU = t_CPU_end - t_CPU_start;

  /* print_vector("w_CPU",h_w,1); */

  fprintf(
          vector_scale_speedup_file,
          "%lf, %lf, %lf,\n",
          t_GPU, t_CPU, t_CPU/t_GPU
          );

  free(active);
  free(h_w);
  free(h_v);

}

void vectorAddOneLoop(FILE* vector_add_speedup_file, const int N, int gridX, int gridY, int threadX, int threadY) {

  int nBytes;
  float *h_w, *h_v, *h_u;

  // Globale Variablen setzen:
  // Anzahl der Inneren Punkte in x- und y-Richtung
  Nx = N;
  Ny = N;

  // Gesamtanzahl der Gitterpunkte
  npts=(Nx+2)*(Ny+2);
  // Aktive Punkte - Array
  active_pts();

  // Speicherbedarf pro Vektor in Byte
  nBytes=npts*sizeof(float);

  // Speicher für Vektoren allozieren
  h_w = (float *) malloc(npts * sizeof(float));
  h_v = (float *) malloc(npts * sizeof(float));
  h_u = (float *) malloc(npts * sizeof(float));

  // auf Null setzen
  memset(h_w, 0, nBytes);
  memset(h_v, 0, nBytes);
  memset(h_u, 0, nBytes);

  // Zufaelliger Vektor
  random_vector(h_v);

  /* print_vector("v",h_v,1); */

  // Device-Speicher allozieren mit cudaMalloc
  float *d_v, *d_w;
  CHECK(cudaMalloc((float**)&d_v, nBytes));
  CHECK(cudaMalloc((float**)&d_w, nBytes));

  // kopieren Host -> Device mit cudaMemcpy
  CHECK(cudaMemcpy(d_v, h_v, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_w, h_w, nBytes, cudaMemcpyHostToDevice));

  dim3 block(gridX, gridY);
  dim3 grid(threadX, threadY);

  double t_GPU_start = seconds();
  laplace2d_GPU <<<block,grid>>> (d_w, d_v, N+2);
  cudaDeviceSynchronize();
  double t_GPU_end = seconds();
  double t_GPU = t_GPU_end - t_GPU_start;

  /* kopieren Device -> Host mit cudaMemcpy */
  CHECK(cudaMemcpy(h_w, d_w, nBytes, cudaMemcpyDeviceToHost));

  /* print_vector("w_GPU",h_w,1); */

  // Device-Speicher freigeben
  CHECK(cudaFree(d_v));
  CHECK(cudaFree(d_w));

  double t_CPU_start = seconds();
  laplace2d_CPU(h_u, h_v);
  double t_CPU_end = seconds();
  double t_CPU = t_CPU_end - t_CPU_start;

  /* print_vector("w_CPU",h_w,1); */

  fprintf(
          vector_add_speedup_file,
          "%lf, %lf, %lf,\n",
          t_GPU, t_CPU, t_CPU/t_GPU
          );

  free(active);
  free(h_w);
  free(h_v);

}

int main(int argc, char **argv)
{
  printf("%s Starting...\n", argv[0]);

  typedef struct grid_parameters {
    int N;
    int gridX;
    int gridY;
    int threadX;
    int threadY;
  } grid_params_t;

  int data_points = 64;

  grid_params_t grid_data[data_points];

  /* read in grid parameters */
  FILE *f = fopen("../scripts/factorizations", "r");
  int i;
  for (i = 0;
       i != data_points &&
         fscanf(f, "%d, %d, %d, %d, %d,\n", &grid_data[i].N, &grid_data[i].gridX, &grid_data[i].gridY, &grid_data[i].threadX, &grid_data[i].threadY) != EOF;
       i++
       );
  fclose(f);

  /* overwrite results */
  FILE *laplace_speedup = fopen("../scripts/laplace_speedup_results", "w");
  FILE *vector_scale_speedup = fopen("../scripts/vector_scale_speedup_results", "w");
  FILE *vector_add_speedup = fopen("../scripts/vector_add_speedup_results", "w");

  fprintf(
          laplace_speedup,
          "t_GPU, t_CPU, t_CPU/t_GPU,\n"
          );
  fprintf(
          vector_scale_speedup,
          "t_GPU, t_CPU, t_CPU/t_GPU,\n"
          );
  fprintf(
          vector_add_speedup,
          "t_GPU, t_CPU, t_CPU/t_GPU,\n"
          );

  int gridX, gridY, threadX, threadY;

  for (i = 0; i < data_points; i++) {

    N = grid_data[i].N;
    gridX = grid_data[i].gridX;
    gridY = grid_data[i].gridY;
    threadX = grid_data[i].threadX;
    threadY = grid_data[i].threadY;

    laplaceOneLoop(
                   laplace_speedup,
                   N,
                   gridX, gridY,
                   threadX, threadY
                   );

    vectorScaleOneLoop(
                       vector_scale_speedup,
                       N,
                       gridX, gridY,
                       threadX, threadY
                       );

    vectorAddOneLoop(
                     vector_add_speedup,
                     N,
                     gridX, gridY,
                     threadX, threadY
                     );

    printf("%d/%d\r", i+1, data_points);
    fflush(stdout);
  }
  printf("\n");

  fclose(laplace_speedup);
  fclose(vector_scale_speedup);
  fclose(vector_add_speedup);

  return (0);

}
