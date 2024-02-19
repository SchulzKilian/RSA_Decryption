#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define NUM_THREADS 8

uint64_t gcd(uint64_t a, uint64_t b) {
    if (b == 0) return a;
    return gcd(b, a % b);
}

uint64_t f(uint64_t x, uint64_t n) {
    return (x * x + 1) % n;
}

uint64_t pollard_rho(uint64_t n) {
    if (n % 2 == 0) return 2;
    uint64_t x = 2, y = 2, d = 1;
    while (d == 1) {
        x = f(x, n);
        y = f(f(y, n), n);
        d = gcd((x > y) ? x - y : y - x, n); // Ensure positive difference
    }
    if (d == n) return 0;
    else return d;
}

typedef struct {
    uint64_t n;
    uint64_t result;
} ThreadData;

void* thread_function(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    data->result = pollard_rho(data->n);
    pthread_exit(NULL);
}

void compute_primes(uint64_t n, uint64_t *result, uint64_t *factor) {
    pthread_t threads[NUM_THREADS];
    ThreadData data[NUM_THREADS];
    int rc;
    long t;

    for (t = 0; t < NUM_THREADS; t++) {
        data[t].n = n;
        rc = pthread_create(&threads[t], NULL, thread_function, (void*)&data[t]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    for (t = 0; t < NUM_THREADS; t++) {
        pthread_join(threads[t], NULL);
        if (data[t].result != 0 && data[t].result != n) {
            *result = data[t].result;
            *factor = (uint64_t)(floor(n / (uint64_t)(*result)));
            break;
        }
    }

    return;
}
