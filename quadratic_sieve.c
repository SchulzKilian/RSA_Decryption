#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <pthread.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cusolverSp.h>
#include <cuda_runtime_api.h>
#include <cusparse_v2.h>


bool is_prime(int n); 
long long mod_exp(long long base, long long exp, long long mod); // Declare mod_exp
void perform_sieving(long long* factor_base, int count, long long n, int*** matrix, int* num_smooth_numbers, int* nnz);
long long* generate_factor_base(long long n, int* count);
cudaError_t transferMatrixToCSRAndToDevice(int** matrix, int numRows, int numCols, int** d_csrRowPtr, int** d_csrColInd, int nnz);






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


#define CHECK_CUDA(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error in file '%s' in line %i: %s.\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
} while (0)

#define CHECK_CUSPARSE(call) do { \
    cusparseStatus_t err = call; \
    if (err != CUSPARSE_STATUS_SUCCESS) { \
        fprintf(stderr, "CUSPARSE error in file '%s' in line %i.\n", __FILE__, __LINE__); \
        exit(EXIT_FAILURE); \
    } \
} while (0)


void solveSparseSystem(int* d_csrRowPtr, int* d_csrColInd, int numRows, int nnz, float* d_b) {
    cusparseHandle_t cusparseH = NULL;
    cusparseSpMatDescr_t matA = NULL;
    cusparseDnVecDescr_t vecX = NULL, vecB = NULL;
    void* dBuffer = NULL;
    float alpha = 1.0;
    float* d_x; 


    CHECK_CUSPARSE(cusparseCreate(&cusparseH));

        
    CHECK_CUDA(cudaMalloc((void**)&d_csrVals, nnz * sizeof(float)));
    float* ones = (float*)malloc(nnz * sizeof(float));
    for (int i = 0; i < nnz; ++i) {
        ones[i] = 1.0;
    }
    CHECK_CUDA(cudaMemcpy(d_csrVals, ones, nnz * sizeof(float), cudaMemcpyHostToDevice));
    free(ones);
    CHECK_CUDA(cudaMalloc((void**)&d_x, numRows * sizeof(float)));


    CHECK_CUSPARSE(cusparseCreateCsr(&matA, numRows, numRows, nnz, d_csrRowPtr, d_csrColInd, d_csrVals,
                                     CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                     CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));

    CHECK_CUSPARSE(cusparseCreateDnVec(&vecB, numRows, d_b, CUDA_R_32F));


    CHECK_CUSPARSE(cusparseCreateDnVec(&vecX, numRows, d_x, CUDA_R_32F));


    size_t bufferSize = 0;
    CHECK_CUSPARSE(cusparseSpSV_bufferSize(cusparseH, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecB, vecX, CUDA_R_32F, CUSPARSE_SPSV_ALG_DEFAULT, spsvDescr, &bufferSize));
    CHECK_CUDA(cudaMalloc(&dBuffer, bufferSize));


    CHECK_CUSPARSE(cusparseSpSV_solve(cusparseH, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecB, vecX, CUDA_R_32F, CUSPARSE_SPSV_ALG_DEFAULT, spsvDescr, dBuffer));


    CHECK_CUDA(cudaFree(dBuffer));
    CHECK_CUSPARSE(cusparseDestroyDnVec(vecX));
    CHECK_CUSPARSE(cusparseDestroyDnVec(vecB));
    CHECK_CUSPARSE(cusparseDestroySpMat(matA));
    CHECK_CUSPARSE(cusparseDestroy(cusparseH));


    float* h_x = (float*)malloc(numRows * sizeof(float));
    cudaMemcpy(h_x, d_x, numRows * sizeof(float), cudaMemcpyDeviceToHost);


    CHECK_CUDA(cudaFree(d_x));

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


void perform_sieving(long long* factor_base, int count, long long n, int*** matrix, int* num_smooth_numbers,int* nnz) {
    *num_smooth_numbers = 0;
    int capacity = n / 100 + 1; // Adjusted initial capacity, ensure it's non-zero

    // Allocate initial space for the matrix
    *matrix = (int**)malloc(capacity * sizeof(int*));
    if (*matrix == NULL) {
        perror("Failed to allocate memory for matrix");
        exit(EXIT_FAILURE);
    }

    // Initialize thread-local storage for all threads
    int max_threads = omp_get_max_threads();
    int*** temp_results = (int***)malloc(max_threads * sizeof(int**)); // Array of pointers to int* arrays
    int* lengths = (int*)calloc(max_threads, sizeof(int)); // Store length of each thread's results
    int* capacities = (int*)malloc(max_threads * sizeof(int)); // Store capacity of each thread's storage

    for (int i = 0; i < max_threads; ++i) {
        capacities[i] = capacity / max_threads + 1; // Initial capacity per thread
        temp_results[i] = (int**)malloc(capacities[i] * sizeof(int*));
        if (temp_results[i] == NULL) {
            perror("Failed to allocate memory for temp_results");
            exit(EXIT_FAILURE);
        }
    }

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        #pragma omp for nowait
        for (long long i = 2; i <= n; ++i) {
            long long num = i;
            int* exponent_vector = (int*)calloc(count, sizeof(int)); // Initialize to 0s
            for (int j = 0; j < count && num > 1; ++j) {
                while (num % factor_base[j] == 0) {
                    exponent_vector[j] = (exponent_vector[j] + 1) % 2; // Calculate exponent modulo 2
                    num /= factor_base[j];
                }
            }
            if (num == 1) { // num is smooth
                if (lengths[id] == capacities[id]) { // Check if the current thread's storage needs expansion
                    capacities[id] *= 2;
                    temp_results[id] = realloc(temp_results[id], capacities[id] * sizeof(int*));
                    if (temp_results[id] == NULL) {
                        perror("Failed to resize temp_results");
                        exit(EXIT_FAILURE);
                    }
                }
                temp_results[id][lengths[id]++] = exponent_vector;
                for (int k = 0; k < count; ++k) {
                    if (exponent_vector[k] != 0) {
                        (*nnz)++;
                    }
                }
            } else {
                free(exponent_vector);
            }
        }
    }

    // Merge results after parallel section
    for (int i = 0; i < max_threads; ++i) {
        for (int j = 0; j < lengths[i]; ++j) {
            if (*num_smooth_numbers == capacity) {
                capacity *= 2; // Resize global matrix if needed
                *matrix = realloc(*matrix, capacity * sizeof(int*));
                if (*matrix == NULL) {
                    perror("Failed to resize matrix");
                    exit(EXIT_FAILURE);
                }
            }
            (*matrix)[(*num_smooth_numbers)++] = temp_results[i][j];
        }
        free(temp_results[i]); // Free each thread's temporary storage
    }
    free(temp_results); // Free the array of pointers
    free(lengths);
    free(capacities);

    // Optionally, resize the matrix to match the exact number of smooth numbers found
    if (*num_smooth_numbers < capacity) {
        *matrix = realloc(*matrix, (*num_smooth_numbers) * sizeof(int*));
        if (*matrix == NULL && *num_smooth_numbers > 0) {
            perror("Failed to shrink memory for matrix");
            // Not exiting here because it's a shrink operation; data is still intact
        }
    }
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
inline cudaError_t checkCuda(cudaError_t result) {
    if (result != cudaSuccess) {
        fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    }
    return result;
}

cudaError_t transferMatrixToCSRAndToDevice(int** denseMatrix, int numRows, int numCols, int** d_csrRowPtr, int** d_csrColInd,  int nnz) {


    // Allocate host memory for CSR format
    int* csrRowPtr = (int*)malloc((numRows + 1) * sizeof(int));
    int* csrColInd = (int*)malloc(nnz * sizeof(int));


    if (!csrRowPtr || !csrColInd) {
        // Handle memory allocation failure
        fprintf(stderr, "Failed to allocate host memory for CSR arrays\n");
        return cudaErrorMemoryAllocation; // Use appropriate CUDA error
    }

    // Fill CSR arrays
    int count = 0;
    csrRowPtr[0] = 0;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            if (denseMatrix[i][j] != 0) {
                csrColInd[count] = j;
                count++;
            }
        }
        csrRowPtr[i + 1] = count;
    }

    printf("CSR Row Pointers:\n");
for (int i = 0; i <= numRows; ++i) {
    printf("%d ", csrRowPtr[i]);
}
printf("\n");

// Print CSR Column Indices
printf("CSR Column Indices:\n");
for (int i = 0; i < nnz; ++i) {
    printf("%d ", csrColInd[i]);
}
printf("\n");
    cudaError_t status;
    status = cudaMalloc((void**)d_csrRowPtr, (numRows + 1) * sizeof(int));
    if (status != cudaSuccess) {
        fprintf(stderr, "Failed to allocate device memory for CSR row pointers\n");
        return status;
    }

    status = cudaMalloc((void**)d_csrColInd, nnz * sizeof(int));
    if (status != cudaSuccess) {
        cudaFree(*d_csrRowPtr);
        fprintf(stderr, "Failed to allocate device memory for CSR column indices\n");
        return status;
    }


    // Copy CSR data from host to device
    status = cudaMemcpy(*d_csrRowPtr, csrRowPtr, (numRows + 1) * sizeof(int), cudaMemcpyHostToDevice);
    if (status != cudaSuccess) {
        fprintf(stderr, "Failed to copy CSR row pointers to device\n");
        goto Error;
    }

    status = cudaMemcpy(*d_csrColInd, csrColInd, nnz * sizeof(int), cudaMemcpyHostToDevice);
    if (status != cudaSuccess) {
        fprintf(stderr, "Failed to copy CSR column indices to device\n");
        goto Error;
    }



    // Free host memory
    free(csrRowPtr);
    free(csrColInd);


    return cudaSuccess;

Error:
    cudaFree(*d_csrRowPtr);
    cudaFree(*d_csrColInd);
    free(csrRowPtr);
    free(csrColInd);
    return status;
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

void factor_primes(long long n) {
    int count;
    long long* factor_base = generate_factor_base(n, &count);
    int** matrix = NULL;
    int num_smooth_numbers = 0;
    int nnz;

    perform_sieving(factor_base, count, n, &matrix, &num_smooth_numbers, &nnz);
    
    int* d_csrRowPtr = NULL;
    int* d_csrColInd = NULL;
    cudaError_t transferStatus = transferMatrixToCSRAndToDevice(matrix, num_smooth_numbers, count, &d_csrRowPtr, &d_csrColInd, nnz);
    if (transferStatus != cudaSuccess) {
        printf("Error transferring matrix to device: %s\n", cudaGetErrorString(transferStatus));
        // Handle error
    } else {
        printf("Matrix transferred to device successfully!\n");
    }
    solveSparseSystem(d_csrRowPtr, d_csrColInd, d_csrVals, num_smooth_numbers, nnz, d_b);

    return 0;
}


int main() {
    long long n;

    // Prompt user for input or set the value of n directly
    printf("Enter a number to factorize: ");
    scanf("%lld", &n);

    // Call the factor_primes function with the given number
    factor_primes(n);

    return 0;
}