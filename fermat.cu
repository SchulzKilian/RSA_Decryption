#include <cuda.h>
#include <stdio.h>
#include <math.h>

__global__ void compute_primes_kernel(long long N, long long* result, long long* factor, int* found) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    long long a = ceil(sqrt((double)N)) + idx;
    long long b2 = a * a - N;
    long long b = llround(sqrt((double)b2));

    if (b * b == b2 && !*found) {
        atomicExch(found, 1);
        atomicExch(result, a - b);
        atomicExch(factor, a + b);
    }
}

void compute_primes(long long N, long long* result, long long* factor) {
    long long *d_result, *d_factor;
    int *d_found;
    int found = 0;
    
    // Allocate memory on device
    cudaMalloc((void**)&d_result, sizeof(long long));
    cudaMalloc((void**)&d_factor, sizeof(long long));
    cudaMalloc((void**)&d_found, sizeof(int));

    // Copy data from host to device
    cudaMemcpy(d_result, result, sizeof(long long), cudaMemcpyHostToDevice);
    cudaMemcpy(d_factor, factor, sizeof(long long), cudaMemcpyHostToDevice);
    cudaMemcpy(d_found, &found, sizeof(int), cudaMemcpyHostToDevice);

    // Launch kernel
    int numBlocks = 1; 
    int threadsPerBlock = 1024; 
    compute_primes_kernel<<<numBlocks, threadsPerBlock>>>(N, d_result, d_factor, d_found);

    // Copy result back to host
    cudaMemcpy(result, d_result, sizeof(long long), cudaMemcpyDeviceToHost);
    cudaMemcpy(factor, d_factor, sizeof(long long), cudaMemcpyDeviceToHost);
    cudaMemcpy(&found, d_found, sizeof(int), cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_result);
    cudaFree(d_factor);
    cudaFree(d_found);
}


