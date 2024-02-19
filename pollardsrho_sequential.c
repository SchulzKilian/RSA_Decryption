#include <stdint.h>

// Iterative version of the GCD function
uint64_t gcd(uint64_t a, uint64_t b) {
    while (b != 0) {
        uint64_t temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

// Pollard's Rho Algorithm for integer factorization remains the same
uint64_t pollards_rho(uint64_t n) {
    uint64_t x = 2, y = 2, d = 1;
    uint64_t c = 1;

    while (d == 1 || d == n || d == 0) {
        x = (x * x + c) % n;
        y = (y * y + c) % n;
        y = (y * y + c) % n;
        d = gcd((x > y) ? x - y : y - x, n);
    }

    return d;
}

// Function to compute primes
void compute_primes(uint64_t n, uint64_t *result, uint64_t *factor) {
    uint64_t factor1, factor2;

    factor1 = pollards_rho(n);
    factor2 = n / factor1;

    *result = factor1;
    *factor = factor2;
    return;
}
