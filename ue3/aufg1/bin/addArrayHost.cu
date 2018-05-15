#include <cuda_runtime.h>
#include <stdio.h>
#include <sys/time.h>

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

void checkResult(float *hostRef, float *gpuRef, const int N)
{
    double epsilon = 1.0E-8;
    bool match = 1;

    for (int i = 0; i < N; i++)
    {
        if (abs(hostRef[i] - gpuRef[i]) > epsilon)
        {
            match = 0;
            printf("Arrays stimmen nicht ueberein!\n");
            printf("host %5.2f gpu %5.2f an der Stelle %d\n", hostRef[i],
                   gpuRef[i], i);
            break;
        }
    }

    if (match) printf("Arrays stimmen ueberein.\n\n");

    return;
}

void initialData(float *ip, int size)
{
    // erzeuge zufaellige Eintraege
    time_t t;
    srand((unsigned) time(&t));

    for (int i = 0; i < size; i++)
    {
        ip[i] = (float)(rand() & 0xFF) / 10.0f;
    }

    return;
}

void addArrayHost(float *A, float *B, float *C, const int N)
{
    for (int idx = 0; idx < N; idx++)
        C[idx] = A[idx] + B[idx];
}

__global__
void addArrayGPU_new(float *A, float *B, float *C) {

  int blockOffset = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y;
  int threadId = threadIdx.x + threadIdx.y * blockDim.x;
  int idx = blockOffset + threadId;

  C[idx] = A[idx] + B[idx];
}

void testGridParamter(const int nElem, FILE* results_file, int blockIdx, int blockIdy, int threadIdx, int threadIdy) {

  // Host-Speicher allozieren mit malloc
  size_t nBytes = nElem * sizeof(float);

  float latenzHDArr, latenzDHArr,
    bandbreiteHDArr, bandbreiteDHArr,
    durchsatzDArr/* , durchsatzHArr */;

  float *h_A, *h_B, *hostRef, *gpuRef;
  h_A     = (float *)malloc(nBytes);
  h_B     = (float *)malloc(nBytes);
  hostRef = (float *)malloc(nBytes);
  gpuRef  = (float *)malloc(nBytes);

  // initialisiere Arrays auf dem Host
  initialData(h_A, nElem);
  initialData(h_B, nElem);

  /* memset(hostRef, 0, nBytes); */
  memset(gpuRef,  0, nBytes);

  // Device-Speicher allozieren mit cudaMalloc
  float *d_A, *d_B, *d_C;
  CHECK(cudaMalloc((float**)&d_A, nBytes));
  CHECK(cudaMalloc((float**)&d_B, nBytes));
  CHECK(cudaMalloc((float**)&d_C, nBytes));

  // Starte Zeitmessung Latenz Host->Device
  double t_HDstart = seconds();

  // kopieren Host -> Device mit cudaMemcpy
  CHECK(cudaMemcpy(d_A, h_A, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_B, h_B, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_C, gpuRef, nBytes, cudaMemcpyHostToDevice));

  // Beende Zeitmessung Latenz Host->Device
  double t_HDend = seconds();

  latenzHDArr= (t_HDend-t_HDstart)*1.e+3;

  // Berechne Bandbreite aus Latenz und Groesse der Arrays
  bandbreiteHDArr = 3*nElem*sizeof(float)/latenzHDArr*1.e-9*1.e+3; // GByte/s


  //Starte Zeitmessung Durchsatz Device
  double t_DurchD_start = seconds();

  dim3 block(blockIdx, blockIdy);
  dim3 grid(threadIdx, threadIdy);
  addArrayGPU_new <<<block, grid>>> (d_A, d_B, d_C);

  cudaDeviceSynchronize();

  //Beende Zeitmessung Durchsatz Device
  double t_DurchD_end = seconds();
  double t_DurchD = t_DurchD_end - t_DurchD_start;

  durchsatzDArr = 3*nElem*1.e-9 /t_DurchD; //Faktor wegen 3 flop

  // Starte Zeitmessung Latenz Device->Host
  double t_DHstart = seconds();

  // kopieren Device -> Host mit cudaMemcpy
  CHECK(cudaMemcpy(gpuRef, d_C, nBytes, cudaMemcpyDeviceToHost));

  // Beende Zeitmessung Latenz Device->Host
  double t_DHend = seconds();

  latenzDHArr = (t_DHend-t_DHstart)*1.e+3;

  // Berechne Bandbreite aus Latenz und Groesse der Arrays
  bandbreiteDHArr = 3*nElem*sizeof(float)/latenzDHArr*1.e-9*1.e+3; // GByte/s
  //printf(" Groesse float %lu \n", sizeof(float));

  /* //Starte Zeitmessung Durchsatz Host */
  /* double t_DurchH_start = seconds(); */

  /* // Addition auf dem Host */
  /* addArrayHost(h_A, h_B, hostRef, nElem); */

  /* //Beende Zeitmessung Durchsatz Host */
  /* double t_DurchH_end = seconds(); */
  /* double t_DurchH = t_DurchH_end - t_DurchH_start; */

  /* durchsatzHArr = nElem*1.e-9 /t_DurchH; */

  /* // verifiziren der Resultate */
  /* checkResult(hostRef, gpuRef, nElem); */

  // Device-Speicher freigeben
  CHECK(cudaFree(d_A));
  CHECK(cudaFree(d_B));
  CHECK(cudaFree(d_C));

  // Host-Speicher freigeben
  free(h_A);
  free(h_B);
  free(hostRef);
  free(gpuRef);

  CHECK(cudaDeviceReset());

  /* write results to file */
  fprintf(
          results_file,
          "%f, %f, %f, %f, %f\n",
          latenzHDArr, latenzDHArr,
          bandbreiteHDArr, bandbreiteDHArr,
          durchsatzDArr/* , durchsatzHArr */
          );
}

int main(int argc, char **argv)
{
  printf("%s Starting...\n", argv[0]);

  typedef struct grid_parameters {
    int gridX;
    int gridY;
    int threadX;
    int threadY;
  } grid_params_t;

  int data_points = 820;
  /* allocate 820 points for the grid data */
  grid_params_t grid_data[data_points];

  /* read in grid parameters */
  FILE *f = fopen("grid_parameters", "r");
  int i;
  for (i = 0;
       i != data_points &&
         fscanf(f, "%d, %d, %d, %d\n", &grid_data[i].gridX, &grid_data[i].gridY, &grid_data[i].threadX, &grid_data[i].threadY) != EOF;
       i++
       );
  fclose(f);

  // Device auswaehlen
  int dev = 0;
  CHECK(cudaSetDevice(dev));

  /* // Groesse der arrays festlegen */
  int nElem = 1024*1024;

  /* overwrite results */
  FILE *g = fopen("grid_parameter_results", "w");

  fprintf(
          g,
          "%s, %s, %s, %s, %s\n",
          "latenzHDArr", "latenzDHArr",
          "bandbreiteHDArr", "bandbreiteDHArr",
          "durchsatzDArr"
          );

  int blockIdx;
  int blockIdy;
  int threadIdx;
  int threadIdy;

  for (i = 0; i < data_points+1; i++) {
    blockIdx = grid_data[i].gridX;
    blockIdy = grid_data[i].gridY;
    threadIdx = grid_data[i].threadX;
    threadIdy = grid_data[i].threadY;

    testGridParamter(nElem, g, blockIdx, blockIdy, threadIdx, threadIdy);
    printf("%d/%d\r", i, data_points);
    fflush(stdout);
  }
  printf("\n");

  /* close results */
  fclose(g);

  return(0);
}
