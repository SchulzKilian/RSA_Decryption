#include <stdint.h>

// Function to compute the greatest common divisor (GCD) of two numbers
uint64_t gcd(uint64_t a, uint64_t b) {
    if (b == 0)
        return a;
    return gcd(b, a % b);
}

// Pollard's Rho Algorithm for integer factorization
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

// Function to compute primess
void compute_primes(uint64_t n, uint64_t *result, uint64_t *factor) {
    uint64_t factor1, factor2;

    factor1 = pollards_rho(n);
    factor2 = n / factor1;

    *result = factor1;
    *factor = factor2;
}
