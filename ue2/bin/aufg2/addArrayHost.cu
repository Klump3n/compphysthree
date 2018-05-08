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

__global__ void addArrayGPU(float *A, float *B, float *C, const int nElem)
{
   if (nElem < 1024) 
   {
   C[threadIdx.x] = A[threadIdx.x] + B[threadIdx.x];  
   }
   
   else 
   {
        int multOfThreads = nElem/1024;
        (int) (nElem % 1024);
        
        for (int k=0; k<multOfThreads; k++)
        {
         if (k*threadIdx.x > nElem) 
         {
             break;
         }
         C[threadIdx.x+1023*k] = A[threadIdx.x+1023*k] + B[threadIdx.x+1023*k];     
        }    
    }
       
}    


int main(int argc, char **argv)
{
    printf("%s Starting...\n", argv[0]);

    // Device auswaehlen
    int dev = 0;
    CHECK(cudaSetDevice(dev));

    // Groesse der arrays festlegen
    int nElem = 2048;
    if (argc>1)
      nElem = atoi(argv[1]);
    printf("Array-Groesse: %d\n", nElem);

    // Host-Speicher allozieren mit malloc
    size_t nBytes = nElem * sizeof(float);

    float *h_A, *h_B, *hostRef, *gpuRef;
    h_A     = (float *)malloc(nBytes);
    h_B     = (float *)malloc(nBytes);
    hostRef = (float *)malloc(nBytes);
    gpuRef  = (float *)malloc(nBytes);

    // initialisiere Arrays auf dem Host
    initialData(h_A, nElem);
    initialData(h_B, nElem);

    memset(hostRef, 0, nBytes);
    memset(gpuRef,  0, nBytes);

    // Device-Speicher allozieren mit cudaMalloc
    float *d_A, *d_B, *d_C;
    CHECK(cudaMalloc((float**)&d_A, nBytes));
    CHECK(cudaMalloc((float**)&d_B, nBytes));
    CHECK(cudaMalloc((float**)&d_C, nBytes));

    // Starte Zeitmessung Latenz Host->Device
    double t_start = seconds();
    
    // kopieren Host -> Device mit cudaMemcpy
    CHECK(cudaMemcpy(d_A, h_A, nBytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_B, h_B, nBytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_C, gpuRef, nBytes, cudaMemcpyHostToDevice));
    
    // Beende Zeitmessung Latenz Host->Device
    double t_end = seconds();
    
    double t_DH = t_end-t_start;
    printf("Kopieren Host -> Device: %f s\n", t_DH);

    // Berechne Durchsatz aus Latenz und Groesse der Arrays
    float durchsatz = 3*nElem*sizeof(float)/t_DH*1.e-9; // GByte/s
    //printf(" Groesse float %lu \n", sizeof(float));
    printf("Durchsatz: %f GB/s\n", durchsatz);
    
    // zu programmieren:
    addArrayGPU<<<1, nElem>>>(d_A, d_B, d_C, nElem);

    // kopieren Host -> Device mit cudaMemcpy
    CHECK(cudaMemcpy(gpuRef, d_C, nBytes, cudaMemcpyDeviceToHost));

    // Addition auf dem Host
    addArrayHost(h_A, h_B, hostRef, nElem);

    // verifiziren der Resultate
    checkResult(hostRef, gpuRef, nElem);

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
    return(0);
}
