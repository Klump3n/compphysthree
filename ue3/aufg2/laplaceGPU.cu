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

  



__global__
void addArrayGPU_new(float *A, float *B, float *C) {

  int blockOffset = (blockIdx.x + blockIdx.y * gridDim.x) * blockDim.x * blockDim.y;
  int threadId = threadIdx.x + threadIdx.y * blockDim.x;
  int idx = blockOffset + threadId;

  C[idx] = A[idx] + B[idx];
}

__global__  
void laplace2d(float *w, float *v, const int N) {
    int i,j;
    int n = N+2;
// i = threadIdx.x, j = threadIdx.y
    i = threadIdx.x + blockIdx.x*blockDim.x;
    j = threadIdx.y + blockIdx.y*blockDim.y;
    printf("(i,j): (%d,%d) \n", i,j);
    
      /* skip the boundary */
      if (
          (i % n == 0) ||
          ((int) (i / n) == 0) ||
          (i % n == n-1) ||
          ((int) (i / n) == n-1)
          ) {
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
  
  
int main(int argc, char **argv)
{
   int N = 1;
   int npts = (N+2)*(N+2);
    
   // Host-Speicher allozieren mit malloc
   size_t nBytes = npts * sizeof(float);
   float *h_v, *h_w;
   h_v = (float *)malloc(nBytes);
   h_w = (float *)malloc(nBytes);
    
   // Device-Speicher allozieren mit cudaMalloc
   float *d_v, *d_w;
   CHECK(cudaMalloc((float**)&d_v, nBytes));
   CHECK(cudaMalloc((float**)&d_w, nBytes));

   // kopieren Host -> Device mit cudaMemcpy
   CHECK(cudaMemcpy(d_v, h_v, nBytes, cudaMemcpyHostToDevice));
   CHECK(cudaMemcpy(d_w, h_w, nBytes, cudaMemcpyHostToDevice));
   
   dim3 block(npts,npts);
   dim3 grid(N,N);
   laplace2d <<<block,grid>>> (d_w, d_v, N);
   cudaDeviceSynchronize();
   
   // kopieren Device -> Host mit cudaMemcpy
   //CHECK(cudaMemcpy(h_w, d_w, nBytes, cudaMemcpyDeviceToHost));  
   
   //print_vector("w",h_w,1);
   
   // Device-Speicher freigeben   
   CHECK(cudaFree(d_v));
   CHECK(cudaFree(d_w));
   
   // Host-Speicher freigeben
   free(h_v);
   free(h_w);
    
}
