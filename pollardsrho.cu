#include <cuda_runtime.h>
#include <stdio.h>

// Utility function: Euclidean algorithm for GCD
__device__ unsigned long long gcd(unsigned long long a, unsigned long long b) {
    while (b != 0) {
        unsigned long long t = b;
        b = a % b;
        a = t;
    }
    return a;
}

// Pollard's Rho CUDA kernel
__global__ void pollardsRhoKernel(unsigned long long N, unsigned long long *d_out, int c) {
    unsigned long long x = 2; // Starting point for x
    unsigned long long y = 2; // Starting point for y
    unsigned long long d = 1; // GCD result initialized to 1

    int threadId = blockIdx.x * blockDim.x + threadIdx.x;
    int c_local = c + threadId; // Adjust c for each thread to vary the function

    while (d == 1) {
        x = (x*x + c_local) % N; // f(x) = (x^2 + c) mod N
        y = (y*y + c_local) % N; // First iteration for y
        y = (y*y + c_local) % N; // Second iteration for y (y moves twice as fast)
        d = gcd(abs(x-y), N);
    }

    if (d != N && d != 1) { // Check if d is a non-trivial factor
        *d_out = d; // Write the non-trivial factor to output
    }
}


void runPollardsRho(unsigned long long N, int c) {
    unsigned long long *d_out; // Device output
    unsigned long long h_out = 0; // Host output

    // Allocat on the device
    cudaMalloc(&d_out, sizeof(unsigned long long));

    // Launch kernel
    int threadsPerBlock = 256;
    int blocks = (N + threadsPerBlock - 1) / threadsPerBlock;

    // Launch the fucking kernel
    pollardsRhoKernel<<<blocks, threadsPerBlock>>>(N, d_out, c);

    // give the result back to the host
    cudaMemcpy(&h_out, d_out, sizeof(unsigned long long), cudaMemcpyDeviceToHost);

    if (h_out > 1) {
        printf("A non-trivial factor of %llu is %llu\n", N, h_out);
    } else {
        printf("No non-trivial factor found for %llu with c = %d\n", N, c);
    }

    // we need device memory
    cudaFree(d_out);
}

int main() {
    unsigned long long N = 589; // Example to factor
    int c = 1; // Constant in f(x) = x^2 + c mod N
    runPollardsRho(N, c);
    return 0;
}
