import random
import intelligent as I
import brute_force as B
import time
import math
import ctypes

# Setup for C function
lib = ctypes.CDLL('./libprimes.so') 
cudalib = ctypes.CDLL('./cudaprimes.so')
lib.compute_primes.argtypes = [ctypes.c_longlong, ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong)]
lib.compute_primes.restype = None
cudalib.compute_primes.argtypes = [ctypes.c_longlong, ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong)]
cudalib.compute_primes.restype = None
result = ctypes.c_longlong()
factor = ctypes.c_longlong()


# Config variables
MIN_RANGE = 10000
MAX_RANGE = 100000
TEST_NUMBERS = 10

def get_prime():
    n = random.choice(primes)
    return n

def test_it():
    cuda = 0
    c_intelligent = 0
    brute_force = 0
    intelligent =  0
    liste = [True,True,True,True]
    for k in range(TEST_NUMBERS):
        b = get_prime()
        c = get_prime()
        if liste[0]:
            start_time= time.time()
            a,d = cudalib.compute_primes(b*c, ctypes.byref(result), ctypes.byref(factor))
            end_time =time.time()
            assert{a,d} == {b,c}
            cuda += end_time - start_time
            print("CUDA  took on average: "+ str(cuda/TEST_NUMBERS)+" seconds to compute")
        if liste[1]:
            start_time= time.time()
            a,d = lib.compute_primes(b*c, ctypes.byref(result), ctypes.byref(factor))
            end_time =time.time()
            assert{a,d} == {b,c}
            c_intelligent += end_time - start_time
            print("C Intelligent  took on average: "+ str(c_intelligent/TEST_NUMBERS)+" seconds to compute")
        if liste[2]:
            start_time= time.time()
            a,d =B.compute_primes(b*c)
            end_time =time.time()
            assert{a,d}== {b,c}
            brute_force += end_time - start_time
            print("Brute Force took on average: "+ str(brute_force/TEST_NUMBERS)+" seconds to compute")
        if liste[3]:
            start_time= time.time()
            a,d = I.compute_primes(b*c)
            end_time =time.time()
            assert{a,d}== {b,c}
            intelligent += end_time - start_time
            print("Intelligent took on average: "+ str(intelligent/TEST_NUMBERS)+" seconds to compute")

    
    
    
    












primes = [160591, 320387, 166609, 652237, 919559, 490151, 528863, 146581, 840187, 924881,
647011, 892733, 960637, 206177, 303011, 354149, 14669, 541579, 600043, 476183,
943079, 877867, 813583, 247613, 209227, 207547, 901253, 569197, 261017, 324791,
697069, 34693, 418351, 734479, 864917, 326119, 603091, 101957, 538397, 990277,
18307, 853187, 593429, 16889, 157901, 675419, 52981, 359641, 44017, 772349,
643231, 587087, 413443, 231551, 497969, 101119, 376841, 586909, 738349, 833377,
723551, 448421, 169489, 416989, 878279, 520381, 900973, 517967, 483649, 85009,
825247, 161573, 268973, 243701, 490591, 870613, 666203, 429413, 118621, 106727,
884243, 960691, 198223, 642361, 139333, 316691, 983441, 904823, 187303, 186437,
770239, 319097, 210127, 887659, 746203, 543503, 86743, 528469, 143609, 593501]




test_it()