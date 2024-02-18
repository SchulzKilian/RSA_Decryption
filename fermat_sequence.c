#include <math.h>
#include <stdio.h>

void compute_primes(long long N, long long* result, long long* factor) {
    long long a = ceil(sqrt(N));


    while (1){
        a++;
        long long b2 = a * a - N;
        long long b = round(sqrt(b2));

        if (b * b == b2) {
            *result = a - b;
            *factor = a + b;
            break; // Found a solution, exit the loop
        }
    }
}
