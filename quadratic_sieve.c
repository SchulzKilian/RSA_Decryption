#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <pthread.h>
#include <math.h>



bool is_prime(int n); 
long long mod_exp(long long base, long long exp, long long mod); // Declare mod_exp
void perform_sieving(long long* factor_base, int count, long long n, int** csrRowPtr, int** csrColInd, int* num_non_zero);
long long* generate_factor_base(long long n, int* count);

void factor_primes(long long n) {
    int count;
    long long* factor_base = generate_factor_base(n, &count);
    int** matrix = NULL;
    int num_smooth_numbers = 0;
    
    // Step 2: sieeeve
    //perform_sieving(factor_base, count, n, &matrix, &num_smooth_numbers);
    
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


void perform_sieving(long long* factor_base, int count, long long n, int** csrRowPtr, int** csrColInd, int* num_non_zero) {
    *num_non_zero = 0;
    *csrRowPtr = (int*)malloc((n + 1) * sizeof(int)); // Allocate memory for CSR row pointer
    *csrColInd = (int*)malloc(n * sizeof(int)); // Allocate memory for CSR column indices

    // Initialize CSR row pointer
    (*csrRowPtr)[0] = 0;

    for (long long i = 2; i <= n; ++i) {
        long long num = i;
        int num_factors = 0;
        for (int j = 0; j < count && num > 1; ++j) {
            while (num % factor_base[j] == 0) {
                num_factors++;
                num /= factor_base[j];
            }
        }
        if (num == 1) { // num is smooth
            for (int j = 0; j < num_factors; ++j) {
                (*csrColInd)[*num_non_zero] = j;
                (*num_non_zero)++;
            }
        }
        (*csrRowPtr)[i] = *num_non_zero; // Update CSR row pointer
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


int countNonZeros(int** hostMatrix, int numRows, int numCols) {
    int count = 0;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            if (hostMatrix[i][j] != 0) {
                ++count;
            }
        }
    }
    return count;
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
        int* csrRowPtr;
    int* csrColInd;
    int num_non_zero;
    // Call the perform_sieving function
    perform_sieving(factor_base, count, n, &csrRowPtr, &csrColInd, &num_non_zero);
    // Print CSR arrays for verification
    printf("CSR Row Pointer:\n");
    for (int i = 0; i <= n; ++i) {
        printf("%d ", csrRowPtr[i]);
    }
    printf("\n\n");
    printf("CSR Column Indices:\n");
    for (int i = 0; i < num_non_zero; ++i) {
        printf("%d ", csrColInd[i]);
    }
    printf("\n");

    // Free allocated memory
    free(csrRowPtr);
    free(csrColInd);
    free(factor_base);
    return 0;
}