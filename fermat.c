#include <math.h>
#include <omp.h>
#include <stdio.h>

void compute_primes(long long N, long long* result, long long* factor) {
    long long a_start = ceil(sqrt(N));
    long long a_end = a_start + 10000; // A predetermined range for demonstration
    int found = 0;

    #pragma omp parallel for
    for (long long a = a_start; a <= a_end; ++a) {
        if (found) continue; // Skip iteration if factors are already found

        long long b2 = a * a - N;
        long long b = round(sqrt(b2));

        if (b * b == b2) {
            #pragma omp critical
            {
                if (!found) { // Ensure only the first solution is taken
                    *result = a - b;
                    *factor = a + b;
                    found = 1;
                }
            }
        }
    }
}
