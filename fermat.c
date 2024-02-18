#include <math.h>
#include <omp.h>
#include <stdio.h>

void compute_primes(long long N, long long* result, long long* factor) {
    long long a_start = ceil(sqrt(N));
    int found = 0;

    // Use a while loop outside the parallel region to control the search range dynamically
    while (!found) {
        long long a_end = a_start + 10; // Adjust the end based on the current start
        
        #pragma omp parallel
        {
            // Each thread checks a portion of the current range
            #pragma omp for
            for (long long a = a_start; a <= a_end; ++a) {
                if (found){
                    continue;
                }
                long long b2 = a * a - N;
                long long b = round(sqrt(b2));
                
                if (b * b == b2) {
                    #pragma omp critical
                    {
                        if (!found) { // Ensure only the first solution is taken
                            *result = a - b;
                            *factor = a + b;
                            found = 1; // Signal to exit the while loop
                        }
                    }
                }
            }
        } // End of parallel region
        
        if (!found) {
            a_start = a_end + 1; // Prepare for the next range if not found
        }
    }
}
