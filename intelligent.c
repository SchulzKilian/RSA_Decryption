#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void compute_primes(long long x, long long *result, long long *factor) {
    if (x % 2 == 0) {
        *result = x / 2;
        *factor = 2;
        return;
    }

    long long length = (long long)sqrt((double)x) + 1;
    for (long long i = 3; i <= length; i += 2) {
        if (x % i == 0) {
            *result = x / i;
            *factor = i;
            return;
        }
    }

    *result = x; // x is prime
    *factor = 1;
}
