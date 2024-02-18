#include <math.h>
#include <omp.h>

void compute_primes(long long N, long long* result, long long* factor) {
    long long a = ceil(sqrt(N));
    long long b2 = a * a - N;
    long long b = round(sqrt(b2));
    int found = 0;

    #pragma omp parallel for private(b2, b) shared(a, found)
    for (long long a = ceil(sqrt(N)); !found; a++) {
        b2 = a * a - N;
        b = round(sqrt(b2));

        if (b * b == b2) {
            #pragma omp critical
            {
                if (!found) { // Ensure only the first solution is taken
                    *factor = a + b;
                    *result = a - b;
                    found = 1; // Signal to other threads that we found the factors
                }
            }
        }
    }
}
