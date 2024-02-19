#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define NUM_THREADS 8

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

// Custom absolute difference function for uint64_t
uint64_t abs_diff(uint64_t a, uint64_t b) {
    return a > b ? a - b : b - a;
}


typedef struct {
    uint64_t n;
    uint64_t seed;
    uint64_t result;
    volatile uint8_t running; // volatile so it can be seen
} ThreadData;

// Now declare pollard_rho_brent with the crrectly defined ThreadData
uint64_t pollard_rho_brent(ThreadData* data);

void* thread_function(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    data->result = pollard_rho_brent(data); // Pass entire data struct
    pthread_exit(NULL);
}

uint64_t pollard_rho_brent(ThreadData* data) {
    uint64_t n = data->n, seed = data->seed;
    if (n % 2 == 0) return 2;
    uint64_t x = seed, y = seed, g = 1, r = 1, q = 1;
    uint64_t ys, m = 128;
    while (g == 1 && data->running) { // Check running flag
        x = y;
        for (uint64_t i = 0; i < r; i++) {
            if (!data->running) return 0; // Early exit if not running
            y = f(y, n);
        }
        uint64_t k = 0;
        while (k < r && g == 1 && data->running) { // Check running flag
            ys = y;
            for (uint64_t i = 0; i < m && i < r - k; i++) {
                if (!data->running) return 0; // Early exit if not running
                y = f(y, n);
                q = modular_mul(q, abs_diff(x, y), n);
            }
            g = gcd(q, n);
            k += m;
        }
        r *= 2;
    }
    if (g == n) return 0; // If g is n, then the factorization failed

    if (g != 1) { // Check if a non-trivial factor was found
        return g;
    }
    
    // If g is 1, try once more with the last ys
    if (!data->running) return 0; // Early exit if not running
    ys = f(ys, n);
    g = gcd(abs_diff(x, ys), n);
    if (g != 1 && g != n) {
        return g;
    }

    return 0; // If no factor was found or if g is n or 1
}


void compute_primes(uint64_t n, uint64_t *result, uint64_t *factor) {
    pthread_t threads[NUM_THREADS];
    ThreadData data[NUM_THREADS];
    int rc;
    long t;

    uint64_t temp_result = 0;
    uint64_t temp_factor = 0;
    int valid_result_found = 0;

    for (t = 0; t < NUM_THREADS; t++) {
        data[t].n = n;
        data[t].seed = 2 + t; // Different seed for each thread
        data[t].running = 1; 
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
            for (int j = 0; j < NUM_THREADS; j++) {
                data[j].running = 0;
            }
            break;
        }
    }

    if (valid_result_found) {
        *result = temp_result;
        *factor = temp_factor;
        return;
    }
}
