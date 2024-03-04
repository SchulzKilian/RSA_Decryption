#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdatomic.h>
#include <time.h>
#include <bits/pthreadtypes.h>

//#define NUM_THREADS 15 // only cause i have 8 cores and for me it works best, different approaches possible
atomic_uint_fast8_t running = 1; 

uint64_t gcd(uint64_t a, uint64_t b) {
    while (b != 0) {
        uint64_t temp = b;
        b = a % b;
        a = temp;
    }

    return a;
}

uint64_t modular_mul(uint64_t a, uint64_t b, uint64_t mod) {
    uint64_t result = 0;
    a %= mod;
    while (b) {
        if (b & 1)
            result = (result + a) % mod;
        a = (2 * a) % mod;
        b >>= 1;
    }
    return result;
}

uint64_t f(uint64_t x, uint64_t n) {
    return (modular_mul(x, x, n) + 1) % n;
}


uint64_t abs_diff(uint64_t a, uint64_t b) {
    return a > b ? a - b : b - a;
}


typedef struct {
    uint64_t n;
    uint64_t seed;
    uint64_t result;

} ThreadData;

// Now declare pollard_rho_brent with the crrectly defined ThreadData
uint64_t pollard_rho_brent(ThreadData* data);

void* thread_function(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    data->result = pollard_rho_brent(data); // Pass entire data struct
    pthread_exit(NULL);
}

uint64_t pollard_rho_brent(ThreadData* data) {
    
    int16_t LOW_GCD_THRESHOLD = 10;
    uint64_t n = data->n, seed = data->seed;
    if (n % 2 == 0) return 2;
    uint64_t x = seed, y = seed, g = 1, r = 1, q = 1;
    uint64_t ys = 0;
    
    uint64_t m = 128;
    while (g == 1 && running) { // Check running flag
        x = y;
        for (uint64_t i = 0; i < r; i++) {
            if (!running) return 0; // Early exit if not running
            y = f(y, n);
        }
        uint64_t k = 0;
        while (k < r && g == 1 && running) { // Check running flag
            ys = y;
            for (uint64_t i = 0; i < m && i < r - k; i++) {
                if (!running) return 0; // Early exit if not running
                y = f(y, n);
                q = modular_mul(q, abs_diff(x, y), n);
            }
            g = gcd(q, n);
            if (g < LOW_GCD_THRESHOLD && g!=1) {
                // This thread found a promising GCD hopefully
                // Adjust parameters to search more aggressively
                printf("It actually got adjusted");
                r *= 10; 
            }

            k += m;
        }
        r *= 2;
    }
    if (g == n) return 0; 

    if (g != 1) { 
        running = 0;
        return g;
    }

    if (!running) return 0;
    ys = f(ys, n);
    g = gcd(abs_diff(x, ys), n);
    if (g != 1 && g != n) {
        return g;
    }

    return 0; // If no factor was found or if g is n or 1
}


void compute_primes(uint64_t n, uint64_t *result, uint64_t *factor, int NUM_THREADS) {
    pthread_t threads[NUM_THREADS];
    ThreadData data[NUM_THREADS];
    int rc;
    long t;


    uint64_t temp_result = 0;
    uint64_t temp_factor = 0;
    int valid_result_found = 0;

    for (t = 0; t < NUM_THREADS; t++) {
        data[t].n = n;
        data[t].seed = 3 + t; // Different seed for each thread

        rc = pthread_create(&threads[t], NULL, thread_function, (void*)&data[t]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    for (t = 0; t < NUM_THREADS; t++) {
        pthread_join(threads[t], NULL);
        if (data[t].result != 0 && data[t].result != n && data[t].result != 1) {
            temp_result = data[t].result;
            temp_factor = n / temp_result;
            valid_result_found = 1;
            break;
        }
    }

    if (valid_result_found) {
        *result = temp_result;
        *factor = temp_factor;
        return;
    }
    
}


int main() {
    uint64_t n = 623546770855596163; // Example number to factorize should be kinda big
    uint64_t result = 0;
    uint64_t factor = 0;
    int thread_counts[] = {1, 2, 4, 8, 15, 20, 25, 32, 64,128}; // Different numbers of threads to test
    int num_configs = sizeof(thread_counts) / sizeof(thread_counts[0]);
    double cpu_time_used[num_configs];

    for (int config = 0; config < num_configs; config++) {
        int NUM_THREADS = thread_counts[config]; 
        double total_time = 0.0;
        int repetitions = 5;

        for (int i = 0; i < repetitions; i++) {
            clock_t start = clock();
            compute_primes(n, &result, &factor, NUM_THREADS); // Adjusted to pass NUM_THREADS
            clock_t end = clock();
            total_time += ((double) (end - start)) / CLOCKS_PER_SEC;
        }

        cpu_time_used[config] = total_time / repetitions;
        printf("Average execution time with %d threads: %f seconds\n", NUM_THREADS, cpu_time_used[config]);
    }

    return 0;
}