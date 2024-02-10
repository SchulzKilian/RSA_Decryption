def compute_primes(x):
    for i in range(2,x):
        result = x/i
        if result == int(result):
            return result, i