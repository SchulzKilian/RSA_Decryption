#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <pthread.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cusparse.h>



bool is_prime(int n); 
long long mod_exp(long long base, long long exp, long long mod); // Declare mod_exp
void perform_sieving(long long* factor_base, int count, long long n, int** csrRowPtr, int* rowPtrSize, int** csrColInd, int* colIndSize);
void printCSRContents(int* csrRowPtr, int rowPtrSize, int* csrColInd, int colIndSize);
long long* generate_factor_base(long long n, int* count);

void factor_primes(long long n) {
    int count;
    long long* factor_base = generate_factor_base(n, &count);
    int** matrix = NULL;
    int num_smooth_numbers = 0;
    
    // Step 2: sieeeve
    perform_sieving(factor_base, count, n, &matrix, &num_smooth_numbers);
    
    // Step 3: and then finish it off
    //solve_linear_algebra(relations, num_relations, count);
    
    // Cleanup
    /*free(factor_base);
    for (int i = 0; i < num_relations; i++) {
        free(relations[i]);
    }
    free(relations);*/
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


void perform_sieving(long long* factor_base, int count, long long n, int** csrRowPtr, int* csrRowPtrSize, int** csrColInd, int* csrColIndSize) {
    int num_smooth_numbers = 0;
    *csrRowPtrSize = 1; // Starting with one element
    *csrColIndSize = 0; // Starting with no elements
    *csrRowPtr = (int*)malloc(*csrRowPtrSize * sizeof(int));
    *csrColInd = NULL; // Will be allocated as needed
    (*csrRowPtr)[0] = 0; // Initial value indicating the start of the first row

    // Initialize thread-local storage for all threads
    int max_threads = omp_get_max_threads();
    int** temp_results = (int**)malloc(max_threads * sizeof(int*));
    int* temp_results_sizes = (int*)calloc(max_threads, sizeof(int));
    int* temp_results_capacities = (int*)malloc(max_threads * sizeof(int));

    for (int i = 0; i < max_threads; ++i) {
        temp_results_capacities[i] = 10; // Arbitrary initial capacity
        temp_results[i] = (int*)malloc(temp_results_capacities[i] * sizeof(int));
    }

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        #pragma omp for nowait
        for (long long i = 2; i <= n; ++i) {
            int* exponent_vector = (int*)calloc(count, sizeof(int)); // Temporarily store indices
            int exponent_vector_size = 0;
            long long num = i;
            for (int j = 0; j < count && num > 1; ++j) {
                while (num % factor_base[j] == 0) {
                    // Instead of pushing back, we manually manage the size and capacity
                    if (exponent_vector_size == temp_results_capacities[id]) {
                        temp_results_capacities[id] *= 2;
                        temp_results[id] = (int*)realloc(temp_results[id], temp_results_capacities[id] * sizeof(int));
                    }
                    temp_results[id][exponent_vector_size++] = j; // Store index
                    num /= factor_base[j];
                }
            }
            if (num == 1) { // num is smooth
                #pragma omp critical
                {
                    for (int k = 0; k < exponent_vector_size; ++k) {
                        if (*csrColIndSize == *csrRowPtrSize - 1) {
                            // Need to expand csrRowPtr and csrColInd
                            *csrRowPtrSize += 1;
                            *csrRowPtr = (int*)realloc(*csrRowPtr, *csrRowPtrSize * sizeof(int));
                            *csrColInd = (int*)realloc(*csrColInd, (*csrColIndSize + exponent_vector_size) * sizeof(int));
                        }
                        (*csrColInd)[(*csrColIndSize)++] = temp_results[id][k];
                    }
                    (*csrRowPtr)[*csrRowPtrSize - 1] = *csrColIndSize;
                    num_smooth_numbers++;
                }
            }
            free(exponent_vector);
        }
    }

    // Clean up
    for (int i = 0; i < max_threads; ++i) {
        free(temp_results[i]);
    }
    free(temp_results);
    free(temp_results_sizes);
    free(temp_results_capacities);
}

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


cudaError_t transferMatrixToDevice(int** hostMatrix, int numRows, int numCols, int** deviceMatrix) {
    int* flatMatrix = NULL;
    size_t size = numRows * numCols * sizeof(int);
    
    // Allocate flat matrix on host
    flatMatrix = (int*)malloc(size);
    if (flatMatrix == NULL) {
        return cudaErrorMemoryAllocation;
    }

    // Flatten the matrix
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            flatMatrix[i * numCols + j] = hostMatrix[i][j];
        }
    }
    return flatMatrix;
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


void print_matrix(int** matrix, int num_smooth_numbers, int count) {
    printf("Matrix of Exponent Vectors Modulo 2:\n");
    for (int i = 0; i < num_smooth_numbers; ++i) {
        for (int j = 0; j < count; ++j) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n"); // New line at the end of each row
    }
}

int main() {
    long long n = 15; // The number to factor
    int count;
    long long* factor_base = generate_factor_base(n, &count);
    int* csrRowPtr = NULL;
    int rowPtrSize = 0;
    int* csrColInd = NULL;
    int colIndSize = 0;

    int num_smooth_numbers = 0;

    perform_sieving(factor_base, count, n, &csrRowPtr, &csrColInd);
    printCSRContents(csrRowPtr,csrColInd)

    int* deviceFlatMatrix;
    copyMatrixToGPU(matrix,num_smooth_numbers,count,&deviceFlatMatrix)
    print_matrix(matrix, num_smooth_numbers, count);

    

    for (int i = 0; i < num_smooth_numbers; ++i) {
        free(matrix[i]); 
    }
    free(matrix); 

    free(factor_base);
    return 0;
}