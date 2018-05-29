#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>

/*
 * This code implements the interleaved Pair approaches to
 * parallel reduction in CUDA. For this example, the sum operation is used.
 */

// Recursive Implementation of Interleaved Pair Approach
int recursiveReduce(int *data, int const size)
{
    // terminate check
    if (size == 1) return data[0];

    // renew the stride
    int const stride = size / 2;

    // in-place reduction
    for (int i = 0; i < stride; i++)
    {
        data[i] += data[i + stride];
    }

    // call recursively
    return recursiveReduce(data, stride);
}

// Kernel: Interleaved Pair Implementation
__global__ void reduceInterleaved (int *g_idata, int *g_odata, unsigned int n)
{
    // set thread ID
    unsigned int tid = threadIdx.x;
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // boundary check
    if(idx >= n) return;

    // in-place reduction in global memory
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2)
    {
        if (tid < stride)
        {
            g_idata[idx] += g_idata[idx + stride];
        }

        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = g_idata[idx];
}

__global__ void reduceUnrolling (int *g_idata, int *g_odata, unsigned int n)
{
    // set thread ID
    unsigned int tid = threadIdx.x;
    unsigned int idx = blockIdx.x * blockDim.x * 2 + threadIdx.x;

    // unroll 2
    if (idx + blockDim.x < n)
    {
        g_idata[idx] += g_idata[idx + blockDim.x];
    }
    __syncthreads();

    // in-place reduction in global memory
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2)
    {
        if (tid < stride)
        {
            g_idata[idx] += g_idata[idx + stride];
        }

        // synchronize within threadblock
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = g_idata[idx];
}

__global__ void qReduceUnrolling (int *g_idata, int *g_odata, unsigned int n, unsigned int q)
{
    // set thread ID
    unsigned int tid = threadIdx.x;
    unsigned int idx = blockIdx.x * blockDim.x * q + threadIdx.x;

    // unroll q
    if (idx + (q-1)*blockDim.x < n)
    {
	for (int k = 1; k < q; k++)
	{
            g_idata[idx] += g_idata[idx + k*blockDim.x ];
	}
    }
    __syncthreads();

    // in-place reduction in global memory
    for (int stride = blockDim.x / 2; stride > 0; stride /= 2)
    {
        if (tid < stride)
        {
            g_idata[idx] += g_idata[idx + stride];
        }

        // synchronize within threadblock
        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = g_idata[idx];
}

__global__
void qReduceUnrollingDouble (double *g_idata, double *g_odata, unsigned int n, unsigned int q)
{
  // set thread ID
  unsigned int tid = threadIdx.x;
  unsigned int idx = blockIdx.x * blockDim.x * q + threadIdx.x;

  // unroll q
  if (idx + (q-1)*blockDim.x < n)
    {
      for (int k = 1; k < q; k++)
        {
          g_idata[idx] += g_idata[idx + k*blockDim.x ];
        }
    }
  __syncthreads();

  // in-place reduction in global memory
  for (int stride = blockDim.x / 2; stride > 0; stride /= 2)
    {
      if (tid < stride)
        {
          g_idata[idx] += g_idata[idx + stride];
        }

      // synchronize within threadblock
      __syncthreads();
    }

  // write result for this block to global mem
  if (tid == 0) g_odata[blockIdx.x] = g_idata[idx];
}

int main(int argc, char **argv)
{

    // set up device
    int dev = 0;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
    printf("%s starting reduction at ", argv[0]);
    printf("device %d: %s ", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));

    bool bResult = false;

    // initialization
    int size = 1 << 24; // total number of elements to reduce
    printf("    with array size %d  ", size);

    // execution configuration
    int blocksize = 512;   // initial block size

    if(argc > 1)
    {
        blocksize = atoi(argv[1]);   // block size from command line argument
    }
    int q;
    if (argc == 3) {
      q = atoi(argv[2]);
    }

    dim3 block (blocksize, 1);
    dim3 grid  ((size + block.x - 1) / block.x, 1);
    printf("grid %d block %d\n", grid.x, block.x);

    // allocate host memory
    size_t bytes = size * sizeof(int);
    int *h_idata = (int *) malloc(bytes);
    int *h_odata = (int *) malloc(grid.x * sizeof(int));
    int *tmp     = (int *) malloc(bytes);

    size_t bytes_d = size * sizeof(double);
    int *h_idata_d = (int *) malloc(bytes_d);
    int *h_odata_d = (int *) malloc(grid.x * sizeof(double));
    int *tmp_d     = (int *) malloc(bytes_d);

    // initialize the array
    int sign=1;
    for (int i = 0; i < size; i++)
    {
        // mask off high 2 bytes to force max number to 255
        h_idata[i] = sign*((int)( rand() & 0xFF ));
        h_idata_d[i] = sign*((double)( rand() & 0xFF ));
        sign*=-1;
    }

    memcpy (tmp, h_idata, bytes);
    memcpy (tmp_d, h_idata_d, bytes_d);

    double iStart, iElaps;
    int gpu_sum = 0;

    // allocate device memory
    int *d_idata = NULL;
    int *d_odata = NULL;
    CHECK(cudaMalloc((void **) &d_idata, bytes));
    CHECK(cudaMalloc((void **) &d_odata, grid.x * sizeof(int)));
    double *d_idata_d = NULL;
    double *d_odata_d = NULL;
    CHECK(cudaMalloc((void **) &d_idata_d, bytes_d));
    CHECK(cudaMalloc((void **) &d_odata_d, grid.x * sizeof(double)));

    // cpu reduction
    iStart = seconds();
    int cpu_sum = recursiveReduce (tmp, size);
    iElaps = seconds() - iStart;
    printf("cpu reduce      elapsed %f sec cpu_sum: %d\n", iElaps, cpu_sum);

    // kernel: reduceInterleaved
    CHECK(cudaMemcpy(d_idata, h_idata, bytes, cudaMemcpyHostToDevice));
    CHECK(cudaDeviceSynchronize());
    iStart = seconds();
    reduceInterleaved<<<grid, block>>>(d_idata, d_odata, size);
    CHECK(cudaDeviceSynchronize());
    iElaps = seconds() - iStart;
    CHECK(cudaGetLastError());
    CHECK(cudaMemcpy(h_odata, d_odata, grid.x * sizeof(int),
                     cudaMemcpyDeviceToHost));
    gpu_sum = 0;

    for (int i = 0; i < grid.x; i++) gpu_sum += h_odata[i];

    printf("gpu Interleaved elapsed %f sec gpu_sum: %d <<<grid %d block "
           "%d>>>\n", iElaps, gpu_sum, grid.x, block.x);

    // kernel: reduceUnrolling
    if (grid.x>1)
    {
       dim3 grid2 ((grid.x + 1)/2,1);
       CHECK(cudaMemcpy(d_idata, h_idata, bytes, cudaMemcpyHostToDevice));
       CHECK(cudaDeviceSynchronize());
       iStart = seconds();
       reduceUnrolling<<<grid2.x, block>>>(d_idata, d_odata, size);
       CHECK(cudaDeviceSynchronize());
       iElaps = seconds() - iStart;
       CHECK(cudaGetLastError());
       CHECK(cudaMemcpy(h_odata, d_odata, grid2.x * sizeof(int),
                        cudaMemcpyDeviceToHost));
       gpu_sum = 0;

       for (int i = 0; i < grid2.x; i++) gpu_sum += h_odata[i];

       printf("gpu Unrolling  elapsed %f sec gpu_sum: %d <<<grid %d block "
              "%d>>>\n", iElaps, gpu_sum, grid2.x, block.x);
    }

    // kernel: q reduceUnrolling
    if (grid.x>1)
    {

      if (!(argc == 3)){
        q = 16; // choose q
      }
       dim3 gridq ((grid.x + 1)/q,1);
       CHECK(cudaMemcpy(d_idata, h_idata, bytes, cudaMemcpyHostToDevice));
       CHECK(cudaDeviceSynchronize());
       iStart = seconds();
       qReduceUnrolling<<<gridq.x, block>>>(d_idata, d_odata, size, q);
       CHECK(cudaDeviceSynchronize());
       iElaps = seconds() - iStart;
       CHECK(cudaGetLastError());
       CHECK(cudaMemcpy(h_odata, d_odata, gridq.x * sizeof(int),
                        cudaMemcpyDeviceToHost));
       gpu_sum = 0;

       for (int i = 0; i < gridq.x; i++) gpu_sum += h_odata[i];




       printf("gpu q Unrolling  elapsed %f sec gpu_sum: %d <<<grid %d block "
              "%d>>>\n", iElaps, gpu_sum, gridq.x, block.x);

       CHECK(cudaMemcpy(d_idata_d, h_idata_d, bytes_d, cudaMemcpyHostToDevice));
       CHECK(cudaDeviceSynchronize());
       iStart = seconds();
       qReduceUnrollingDouble<<<gridq.x, block>>>(d_idata_d, d_odata_d, size, q);
       CHECK(cudaDeviceSynchronize());
       iElaps = seconds() - iStart;
       CHECK(cudaGetLastError());
       CHECK(cudaMemcpy(h_odata_d, d_odata_d, gridq.x * sizeof(double),
                        cudaMemcpyDeviceToHost));
       gpu_sum = 0;

       for (int i = 0; i < gridq.x; i++) gpu_sum += h_odata[i];
    }

    // free host memory
    free(h_idata);
    free(h_odata);
    free(h_idata_d);
    free(h_odata_d);

    // free device memory
    CHECK(cudaFree(d_idata));
    CHECK(cudaFree(d_odata));
    CHECK(cudaFree(d_idata_d));
    CHECK(cudaFree(d_odata_d));

    // reset device
    CHECK(cudaDeviceReset());

    // check the results
    bResult = (gpu_sum == cpu_sum);

    if(!bResult) printf("Test failed!\n");

    return EXIT_SUCCESS;
}
