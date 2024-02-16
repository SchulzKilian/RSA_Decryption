#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <pthread.h>
#include <math.h>


bool is_prime(int n); 
long long mod_exp(long long base, long long exp, long long mod); // Declare mod_exp
void perform_sieving(long long* factor_base, int count, long long n, int*** matrix, int* num_smooth_numbers);
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


void perform_sieving(long long* factor_base, int count, long long n, int*** matrix, int* num_smooth_numbers) {
    *num_smooth_numbers = 0;
    int capacity = n/100 + 1; // Adjusted initial capacity

    // Allocate initial space for the matrix
    *matrix = (int**)malloc(capacity * sizeof(int*));
    if (*matrix == NULL) {
        perror("Failed to allocate memory for matrix");
        exit(EXIT_FAILURE);
    }

    // Temporary storage for thread-local results
    int** temp_matrix = malloc(capacity * sizeof(int*));
    if (temp_matrix == NULL) {
        perror("Failed to allocate memory for temp_matrix");
        exit(EXIT_FAILURE);
    }
    int temp_num_smooth = 0;

    #pragma omp parallel
    {
        int** local_matrix = NULL;
        int local_num_smooth = 0;
        int local_capacity = capacity;

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
                if (local_num_smooth == local_capacity) {
                    // Resize local_matrix if needed, similar to realloc logic
                    local_capacity *= 2;
                    local_matrix = realloc(local_matrix, local_capacity * sizeof(int*));
                    if (local_matrix == NULL) {
                        perror("Failed to resize local_matrix");
                        exit(EXIT_FAILURE);
                    }
                }
                local_matrix[local_num_smooth++] = exponent_vector;
            } else {
                free(exponent_vector);
            }
        }

        #pragma omp critical
        {
            for (int i = 0; i < local_num_smooth; ++i) {
                if (*num_smooth_numbers >= capacity) {
                    // Resize global matrix if needed, similar to realloc logic
                    capacity *= 2;
                    *matrix = realloc(*matrix, capacity * sizeof(int*));
                    if (*matrix == NULL) {
                        perror("Failed to resize matrix");
                        exit(EXIT_FAILURE);
                    }
                }
                (*matrix)[(*num_smooth_numbers)++] = local_matrix[i];
            }
            free(local_matrix);
        }
    }

    free(temp_matrix);

    // Optional: Resize the matrix to match the exact number of smooth numbers found
    *matrix = realloc(*matrix, (*num_smooth_numbers) * sizeof(int*));
    if (*matrix == NULL && *num_smooth_numbers > 0) {
        perror("Failed to shrink memory for matrix");
        exit(EXIT_FAILURE);
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
    long long n = 221; // The number to factor
    int count;
    long long* factor_base = generate_factor_base(n, &count);
    int** matrix = NULL;
    int num_smooth_numbers = 0;

    perform_sieving(factor_base, count, n, &matrix, &num_smooth_numbers);


    print_matrix(matrix, num_smooth_numbers, count);

    

    for (int i = 0; i < num_smooth_numbers; ++i) {
        free(matrix[i]); 
    }
    free(matrix); 

    free(factor_base);
    return 0;
}