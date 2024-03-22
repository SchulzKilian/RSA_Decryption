#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <stdatomic.h>
#include <time.h>
#include <bits/pthreadtypes.h>

//atomic_uint_fast8_t running = 1; 
atomic_bool running = ATOMIC_VAR_INIT(1);

void gcd(mpz_t result, const mpz_t a, const mpz_t b) {
    mpz_t tmp_a, tmp_b;
    mpz_init_set(tmp_a, a);
    mpz_init_set(tmp_b, b);
    while (mpz_cmp_ui(tmp_b, 0) != 0) {
        mpz_t temp;
        mpz_init_set(temp, tmp_b);
        mpz_mod(tmp_b, tmp_a, tmp_b);
        mpz_set(tmp_a, temp);
        mpz_clear(temp);
    }
    mpz_set(result, tmp_a);
    mpz_clear(tmp_a);
    mpz_clear(tmp_b);
}

void modular_mul(mpz_t result, const mpz_t a, const mpz_t b, const mpz_t mod) {
    mpz_t tmp_a, tmp_b;
    mpz_init_set(tmp_a, a);
    mpz_init_set(tmp_b, b);
    mpz_mod(tmp_a, tmp_a, mod);
    mpz_set_ui(result, 0);
    while (mpz_cmp_ui(tmp_b, 0) > 0) {
        if (mpz_odd_p(tmp_b))
            mpz_add(result, result, tmp_a);
        mpz_mod(result, result, mod);
        mpz_mul_ui(tmp_a, tmp_a, 2);
        mpz_mod(tmp_a, tmp_a, mod);
        mpz_fdiv_q_2exp(tmp_b, tmp_b, 1); // b >>= 1
    }
    mpz_mod(result, result, mod);
    mpz_clear(tmp_a);
    mpz_clear(tmp_b);
}

void f(mpz_t result, const mpz_t x, const mpz_t n) {
    mpz_t tmp;
    mpz_init(tmp);
    modular_mul(tmp, x, x, n);
    mpz_add_ui(tmp, tmp, 1);
    mpz_mod(result, tmp, n);
    mpz_clear(tmp);
}

void abs_diff(mpz_t result, const mpz_t a, const mpz_t b) {
    if (mpz_cmp(a, b) > 0) {
        mpz_sub(result, a, b);
    } else {
        mpz_sub(result, b, a);
    }
}

typedef struct {
    mpz_t n;
    unsigned long int seed;
    mpz_t result;
} ThreadData;

void pollard_rho_brent(ThreadData* data);

void* thread_function(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    pollard_rho_brent(data); // Pass entire data struct
    pthread_exit(NULL);
}

void pollard_rho_brent(ThreadData* data) {
    mpz_t one, two;
    mpz_init_set_ui(one, 1);
    mpz_init_set_ui(two, 2);
    
    mpz_t n, x, y, g, r, q, ys, m, i, k, tmp_diff;
    mpz_inits(n, x, y, g, r, q, ys, m, i, k, tmp_diff, NULL);
    mpz_set(n, data->n);
    
    if (mpz_even_p(n)) {
        mpz_set(data->result, two);
        mpz_clears(one, two, n, x, y, g, r, q, ys, m, i, k, tmp_diff, NULL);
        return;
    }

    mpz_set_ui(x, data->seed);
    mpz_set_ui(y, data->seed);
    mpz_set_ui(g, 1);
    mpz_set_ui(r, 1);
    mpz_set_ui(q, 1);
    mpz_set_ui(m, 128);
    int16_t LOW_GCD_THRESHOLD = 10;

    while (!mpz_cmp_ui(g, 1) && atomic_load(&running)) {
        mpz_set(x, y);
        for (mpz_set_ui(i, 0); mpz_cmp(i, r) < 0; mpz_add_ui(i, i, 1)) {
            f(y, y, n);
        }
        mpz_set_ui(k, 0);
        while (mpz_cmp(k, r) < 0 && !mpz_cmp_ui(g, 1) && atomic_load(&running)) {
            mpz_set(ys, y);
            for (mpz_set_ui(i, 0); mpz_cmp(i, m) < 0 && mpz_cmp(i, r) < 0; mpz_add_ui(i, i, 1)) {
                f(y, y, n);
                abs_diff(tmp_diff, x, y);
                modular_mul(q, q, tmp_diff, n);
            }
            gcd(g, q, n);
            mpz_add(k, k, m);
        }
        mpz_mul_ui(r, r, 2);
    }
    if (mpz_cmp(g, n) == 0 || !atomic_load(&running)) {
        mpz_set_ui(data->result, 0);
    } else if (mpz_cmp_ui(g, 1)) {
        atomic_store(&running, 0);
        mpz_set(data->result, g);
    }

    mpz_clears(one, two, n, x, y, g, r, q, ys, m, i, k, tmp_diff, NULL);
}

void compute_primes(mpz_t n, mpz_t *result, mpz_t *factor, int NUM_THREADS) {
    pthread_t threads[NUM_THREADS];
    ThreadData data[NUM_THREADS];
    int rc;

    // Initialize results and factors with 0 values
    mpz_init(*result);
    mpz_init(*factor);
    int valid_result_found = 0;
    running = 1;
    int warp = 0;
    while (valid_result_found == 0){
    for (int t = 0; t < NUM_THREADS; t++) {
        mpz_init(data[t].n);
        mpz_init(data[t].result);
        mpz_set(data[t].n, n);
        data[t].seed = 2 + t+NUM_THREADS*warp; // Different seed for each thread

        rc = pthread_create(&threads[t], NULL, thread_function, (void *)&data[t]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    
    for (int t = 0; t < NUM_THREADS; t++) {
        pthread_join(threads[t], NULL);
        if (mpz_cmp_ui(data[t].result, 0) != 0 && mpz_cmp(data[t].result, n) != 0) {
            mpz_set(*result, data[t].result);
            mpz_divexact(*factor, n, *result); // *factor = n / *result;
            valid_result_found = 1;
            break;
        }

    }


    // Clean up
    for (int t = 0; t < NUM_THREADS; t++) {
        mpz_clear(data[t].n);
        mpz_clear(data[t].result);
    }
    warp++;
    }
    if (!valid_result_found) {
        mpz_set_ui(*result, 0);
        mpz_set_ui(*factor, 0);
    }
}


int main() {
    char* semiprimes[] = {
        "281903819787774937209232578689",
        "364484869645843621622832704153",
        "345469451591760174730230603691",
        "463398692478296989516736649887",
        "332227475082334623786092387273",
        "31626851437447153694819292167",
        "777017557862414631890710147753",
        "147496354805768937296258608211",
        "538232757297255111133394171213",
        "137486608107202789957046404849"
    };
    int num_semiprimes = sizeof(semiprimes) / sizeof(semiprimes[0]);
    char* medium_primes[] = {
        "2080615832610841335480683",
        "1438659423804570836893217",
        "1111321301231058844462511",
        "744777372083245989163421",
        "6053447060883065240701121",
        "3584405191726040249488253",
        "2682848904165256275971861",
        "4225513033754881584728239",
        "4605131422979549619984269",
        "5278982110844980463386933"
    };
    int num_medium_primes = sizeof(medium_primes) / sizeof(medium_primes[0]);
        char* small_semiprimes[] = {
        "2398115620361697",
        "2821095062496139",
        "3259994183448281",
        "3714900574054983",
        "4185934403316433",
        "4673326814780727",
        "5177308227255299",
        "5698119380409281",
        "6235999965097863",
        "6791200376458947"
    };

    int num_small_semiprimes = sizeof(small_semiprimes) / sizeof(small_semiprimes[0]);
    int threadnumbers[7] = {4,2,1};
    int num_threads = sizeof(threadnumbers) / sizeof(threadnumbers[0]);

    for (int i = 0; i < num_threads; i++) {
        int NUM_THREADS = threadnumbers[i];
        double total_time_taken = 0.0;

        for (int j = 0; j < 1; j++) {
            mpz_t n, result, factor;
            mpz_init_set_str(n, medium_primes[j], 10); // Initialize n with a large number
            mpz_init(result);
            mpz_init(factor);
            atomic_store(&running, 1);

            clock_t start = clock();
            compute_primes(n, &result, &factor, NUM_THREADS);
            clock_t end = clock();

            if (mpz_cmp_ui(result, 0) != 0) { 
                gmp_printf("Prime factor of %Zd: %Zd and %Zd\n", n, result, factor);
            } else {
                gmp_printf("No prime factor found for %Zd.\n", n);
            }

            double time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
            total_time_taken += time_taken;

            mpz_clear(n);
            mpz_clear(result);
            mpz_clear(factor);
        }
        
        double average_time_taken = total_time_taken / num_medium_primes;
        printf("Average execution time for %d threads: %.2f seconds\n", NUM_THREADS, average_time_taken);
    }

    return 0;
}