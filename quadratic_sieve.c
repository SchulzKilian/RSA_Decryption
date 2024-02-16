#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <pthread.h>
#include <math.h>
//#include <cuda_runtime.h>
bool is_prime(int n); // Assuming this is defined correctly elsewhere
long long mod_exp(long long base, long long exp, long long mod); // Declare mod_exp

long long* generate_factor_base(long long n, int* count);

void factor_primes(long long n) {
    int count; // Number of primes in the factor base
    int num_relations; // Number of relations found during sieving
    
    // Step 1: get the factor thind
    long long* factor_base = generate_factor_base(n, &count);
    
    // Step 2: sieeeve
    //int* relations = perform_sieving(factor_base, count, n, &num_relations);
    
    // Step 3: and then finish it off
    //solve_linear_algebra(relations, num_relations, count);
    
    // Cleanup
    //free(factor_base);
    //for (int i = 0; i < num_relations; i++) {
    //    free(relations[i]);
    //}
    //free(relations);
}




// very inefficient will change soon to eristo something prime generation
bool is_prime(int n) {
    if (n <= 1) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;
    for (int i = 3; i <= sqrt(n); i += 2) {
        if (n % i == 0) return false;
    }
    return true;
}

long long* generate_factor_base(long long n, int* count) {
    int size = (int)exp(sqrt(log(n) * log(log(n))));   // a common heuristic for how many prime numbers to check
    
    long long* factor_base = malloc(size * sizeof(long long));
    if (!factor_base) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    int actual_count = 0;
    
    #pragma omp parallel
    {
        int local_count = 0;
        long long* local_factor_base = malloc(size * sizeof(long long));
        
        #pragma omp for nowait
        for (int i = 2; i <= size; i++) {
            if (is_prime(i) && mod_exp(n, (i - 1) / 2, i) == 1) { // Assuming the thing  is defined elsewhere
                local_factor_base[local_count++] = i;
            }
        }
        
        #pragma omp critical
        {
            for (int j = 0; j < local_count; j++) {
                factor_base[actual_count++] = local_factor_base[j];
            }
        }
        
        free(local_factor_base);
    }
    
    *count = actual_count;
    return factor_base;
}


/*void perform_sieving(long long* factor_base, int factor_base_size, int range_start, int range_size) {
    long long *d_factor_base;
    int *d_relations;

    // Allocate GPU memory
    cudaMalloc(&d_factor_base, factor_base_size * sizeof(long long));
    cudaMalloc(&d_relations, range_size * factor_base_size * sizeof(int));

    // Copy factor base to GPU
    cudaMemcpy(d_factor_base, factor_base, factor_base_size * sizeof(long long), cudaMemcpyHostToDevice);

    // Kernel launch parameters
    int threadsPerBlock = 256;
    int blocksPerGrid = (range_size + threadsPerBlock - 1) / threadsPerBlock;

    // Launch the kernel
    sieve_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_factor_base, factor_base_size, d_relations, range_start, range_size);

    // Allocate host memory for relations and copy back from GPU
    int* relations = (int*)malloc(range_size * factor_base_size * sizeof(int));
    cudaMemcpy(relations, d_relations, range_size * factor_base_size * sizeof(int), cudaMemcpyDeviceToHost);

    // Process the relations here...
    // Example: Print out relations for each smooth number
    for (int i = 0; i < range_size; ++i) {
        std::cout << "Relations for number " << (range_start + i) << ": ";
        for (int j = 0; j < factor_base_size; ++j) {
            std::cout << relations[i * factor_base_size + j] << " ";
        }
        std::cout << "\n";
    }

    // Cleanup
    cudaFree(d_factor_base);
    cudaFree(d_relations);
    free(relations);
}


__global__ void sieve_kernel(long long *factor_base, int factor_base_size, int *relations, int range_start, int range_size) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    for (int i = index; i < range_size; i += stride) {
        int num = range_start + i;
        for (int j = 0; j < factor_base_size; j++) {
            if (num % factor_base[j] == 0) {
                // Mark this prime as dividing the number
                relations[i * factor_base_size + j] = 1;
            } else {
                relations[i * factor_base_size + j] = 0;
            }
        }
    }
}

*/

long long mod_exp(long long base, long long exp, long long mod) {
    long long result = 1;
    base = base % mod;
    while (exp > 0) {
        // If exp is odd, multiply the base with result
        if (exp % 2 == 1)
            result = (result * base) % mod;
        // Now exp must be even, square the base
        exp = exp >> 1; // Equivalent to exp = exp / 2
        base = (base * base) % mod;
    }
    return result;
}

int main() {
    long long n = 84923; // The number to factor
    int count;
    long long* factor_base = generate_factor_base(n, &count);

    printf("Factor base for %lld:\n", n);
    for (int i = 0; i < count; i++) {
        printf("%lld ", factor_base[i]);
    }
    printf("\nCount: %d\n", count);

    free(factor_base);
    return 0;
}