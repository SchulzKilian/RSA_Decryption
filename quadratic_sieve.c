#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <pthread.h>
#include <math.h>


bool is_prime(int n); 
long long mod_exp(long long base, long long exp, long long mod); // Declare mod_exp
void perform_sieving(long long* factor_base, int count, long long n, int* num_relations);
long long* generate_factor_base(long long n, int* count);

void factor_primes(long long n) {
    int count; // Number of primes in the factor base
    int num_relations; // Number of relations found during sieving
    
    // Step 1: get the factor thind
    long long* factor_base = generate_factor_base(n, &count);
    
    // Step 2: sieeeve
    perform_sieving(factor_base, count, n, &num_relations);
    
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


void perform_sieving(long long* factor_base, int count, long long n, int* num_relations) {
    // Initialize num_relations to 0
    *num_relations = 0;

    // Allocate an array thing where i can check arew they smooth
    bool* isCandidate = (bool*)calloc(n + 1, sizeof(bool)); // Initialized to true, easier

    // caan i factor it?
    #pragma omp parallel for schedule(dynamic)
    for (long long i = 2; i <= n; ++i) {
        long long num = i;
        for (int j = 0; j < count && num > 1; ++j) {
            long long prime = factor_base[j];
            while (num % prime == 0) {
                num /= prime; // Divide by prime as much as possible
            }
        }
        // If num is reduced to 1, all its factors are in the factor base
        if (num == 1) {
            #pragma omp atomic
            (*num_relations)++;
            #pragma omp critical
            {
                printf("%lld ", i); // Print smooth numbers
            }
        }
    }

    free(isCandidate);
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

int main() {
    long long n = 30; // The number to factor
    int count;
    long long* factor_base = generate_factor_base(n, &count);

    printf("Factor base for %lld:\n", n);
    for (int i = 0; i < count; i++) {
        printf("%lld ", factor_base[i]);
    }
    printf("\nCount of primes in factor base: %d\n", count);

    // Additions for perform_sieving verification
    int num_relations;
    perform_sieving(factor_base, count, n, &num_relations);
    printf("Number of smooth numbers up to %lld: %d\n", n, num_relations);


    if (num_relations == 7) {
        printf("The sieving process found the right amount ofsmooth numbers, which suggests it is working.\n");
    } else {
        printf("Not the right smooth numbers found. There may be an issue with the sieving process or the factor base.\n");
    }

    free(factor_base);
    return 0;
}