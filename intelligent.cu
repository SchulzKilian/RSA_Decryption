#include <stdio.h>
#include <cuda_runtime.h>

__global__ void compute_primes_kernel(long long x, long long* result) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    long long number_to_check = idx + 2; // Start checking from 2

    if (number_to_check * number_to_check <= x) {
        if (x % number_to_check == 0) {
            result[0] = x / number_to_check;
            result[1] = number_to_check;
        }
    }
}

extern "C" void compute_primes(long long x, long long* host_result) {
    long long* device_result;
    cudaMalloc((void**)&device_result, 2 * sizeof(long long));
    cudaMemset(device_result, 0, 2 * sizeof(long long));

    int threadsPerBlock = 256;
    int blocksPerGrid = (x + threadsPerBlock - 1) / threadsPerBlock;
    compute_primes_kernel<<<blocksPerGrid, threadsPerBlock>>>(x, device_result);

    cudaMemcpy(host_result, device_result, 2 * sizeof(long long), cudaMemcpyDeviceToHost);
    cudaFree(device_result);
}

int main() {
    long long x = 22742734291LL * 52711LL;
    long long host_result[2] = {0};

    compute_primes(x, host_result);

    if (host_result[0] != 0) {
        printf("Factor: %lld, Result: %lld\n", host_result[1], host_result[0]);
    } else {
        printf("No factors found.\n");
    }

    return 0;
}
