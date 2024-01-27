#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void compute_primes(long long x, long long *result, long long *factor) {
    long long length = (long long)sqrt((double)x) + 1;
    int *l = (int*)malloc(length * sizeof(int));

    if (l == NULL) {
        *result = -1; // Memory allocation error
        return;
    }

    for (long long i = 0; i < length; i++) {
        l[i] = 0;
    }

    for (long long ind = 0; ind < length; ind++) {
        if (l[ind] || ind == 0)
            continue;

        long long numb = ind + 1;
        double temp_result = (double)x / numb;

        if (temp_result == (long long)temp_result) {
            *result = (long long)temp_result;
            *factor = numb;
            free(l);
            return;
        } else {
            l[ind] = 1;
            long long k = numb;
            while (k * numb <= length - 1 && k < 10) {
                l[k * numb - 1] = 1;
                k += numb;
            }
        }
    }

    *result = 0; // No factors found
    free(l);
}