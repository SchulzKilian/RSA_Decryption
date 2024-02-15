This is an attempt at a parallelized and fast way to find the two factors of a large number that has two prime factors, as needed for the RSA algorithm.
The goal is for the parallelization to be three times as fast as the brute force in normal c.

Different algorithmic attempts:

- brute force (quite fast and can be used for data parallelization)
- brute force only trying with numbers where the last bit is not zero (almost double as fast)
- initializing a boolean array and then cancelling out multiples of checked numbers (slower than brute force)
- brute force but checking if number is prime instead of divisor (terrible idea idk why i tried)
- fermats factorization (only works if numbers are close to eachother)
- pollards rho (trying to implement to test)
- quadratic sieve (efficient for """smaller""" numbers, sufficient for what i am trying to reach)
- general number field sieve (didnt try yet because it seems super complicated, will be last resort)