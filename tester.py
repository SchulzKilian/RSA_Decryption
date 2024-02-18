import random
import intelligent as I
import brute_force as B

import time
import math
import ctypes

# Config variables
MIN_RANGE = 10000
MAX_RANGE = 100000
TEST_NUMBERS = 5

liste = [
         False, # cuda normal
         False, # c normal
         False, # brute force pythong
         False, # sieve python
         False, # fermat openmp
         True,# Quadratic Sieve
         True # fermat sequence
         ]







# Setup for C functions
if liste[1]:
    lib = ctypes.CDLL('./libprimes.so') 
    lib.compute_primes.argtypes = [ctypes.c_longlong, ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong)]
    lib.compute_primes.restype = None
    result = ctypes.c_longlong()
    factor = ctypes.c_longlong()

if liste[0]:
    cudalib = ctypes.CDLL('./cudaprimes.so')
    cudalib.compute_primes.argtypes = [ctypes.c_longlong, ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong)]
    cudalib.compute_primes.restype = None
    result = ctypes.c_longlong()
    factor = ctypes.c_longlong()
if liste[4]:
    fermatlib = ctypes.CDLL('./libfermat.so') 
    fermatlib.compute_primes.argtypes = [ctypes.c_longlong, ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong)]
    fermatlib.compute_primes.restype = None
    result = ctypes.c_longlong()
    factor = ctypes.c_longlong()

if liste[5]:
    quadlib = ctypes.CDLL('./libfermat.so') 
    quadlib.compute_primes.argtypes = [ctypes.c_longlong, ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong)]
    quadlib.compute_primes.restype = None
    result = ctypes.c_longlong()
    factor = ctypes.c_longlong()

if liste[6]:
    fermatser = ctypes.CDLL('./fermat_sequence.so') 
    fermatser.compute_primes.argtypes = [ctypes.c_longlong, ctypes.POINTER(ctypes.c_longlong), ctypes.POINTER(ctypes.c_longlong)]
    fermatser.compute_primes.restype = None
    result = ctypes.c_longlong()
    factor = ctypes.c_longlong()   



def get_prime():
    n = random.choice(primes)
    return n

def get_medium_small_prime():
    n = random.choice(medium_small_primes)
    return n

def get_medium_prime():
    n = random.choice(medium_primes)
    return n

def get_big_prime():
    n = random.choice(big_primes)
    return n

def test_it():
    cuda = 0
    c_intelligent = 0
    brute_force = 0
    fermat = 0
    intelligent =  0
    fermats = 0
    q_sieve= 0
    for k in range(TEST_NUMBERS):
        b = get_medium_prime()
        c = get_medium_prime()
        if liste[0]:
            start_time= time.time()
            a,d = cudalib.compute_primes(b*c, ctypes.byref(result), ctypes.byref(factor))
            end_time =time.time()
            assert{a,d} == {b,c}
            cuda += end_time - start_time
            if k == TEST_NUMBERS-1:
                print("CUDA  took on average: "+ str(cuda/TEST_NUMBERS)+" seconds to compute")
        if liste[1]:
            start_time= time.time()
            lib.compute_primes(b*c, ctypes.byref(result), ctypes.byref(factor))
            a,d = int(result.value), int(factor.value)
            end_time =time.time()
            try:
                assert{a,d} == {b,c}
            except:
                print(a,d)
                print(b,c)
            c_intelligent += end_time - start_time
            if k == TEST_NUMBERS-1:
                print("C Intelligent  took on average: "+ str(c_intelligent/TEST_NUMBERS)+" seconds to compute")
        if liste[2]:
            start_time= time.time()
            a,d =B.compute_primes(b*c)
            end_time =time.time()
            assert{a,d}== {b,c}
            brute_force += end_time - start_time
            if k == TEST_NUMBERS-1:
                print("Brute Force took on average: "+ str(brute_force/TEST_NUMBERS)+" seconds to compute")
        if liste[3]:
            start_time= time.time()
            a,d = I.compute_primes(b*c)
            end_time =time.time()
            assert{a,d}== {b,c}
            intelligent += end_time - start_time
            if k == TEST_NUMBERS-1:
                print("Intelligent took on average: "+ str(intelligent/TEST_NUMBERS)+" seconds to compute")
        if liste[4]:
            start_time= time.time()
            fermatlib.compute_primes(b*c, ctypes.byref(result), ctypes.byref(factor))
            a,d = int(result.value), int(factor.value)
            end_time =time.time()
            assert{a,d}== {b,c}
            fermat += end_time - start_time
            if k == TEST_NUMBERS-1:
                print("Fermat took on average: "+ str(fermat/TEST_NUMBERS)+" seconds to compute")
        if liste[5]:
            start_time= time.time()
            lib.factor_primes(b*c)
            a,d = int(result.value), int(factor.value)
            end_time =time.time()
            q_sieve += end_time - start_time
            if k == TEST_NUMBERS-1:
                print("The quadratic sieving took on average: "+ str(q_sieve/TEST_NUMBERS)+" seconds to compute")
        if liste[6]:
            start_time= time.time()
            fermatser.compute_primes(b*c, ctypes.byref(result), ctypes.byref(factor))
            a,d = int(result.value), int(factor.value)
            end_time =time.time()
            try:
                assert{a,d} == {b,c}
            except:
                print(a,d)
                print(b,c)
            fermats += end_time - start_time
            if k == TEST_NUMBERS-1:
                print("Fermat in sequence took on average: "+ str(fermats/TEST_NUMBERS)+" seconds to compute")

    
    
    
    




 







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


big_primes = [
    316691868537677,
    212588137232143,
    931248599717287,
    462481314583859,
    427965809376023,
    357900188652949,
    469738788499069,
    968377619835011,
    501361019742467,
    874776702088369,
    549297986383687,
    916679942097833,
    363173953821421,
    582527240138353,
    663572060523377,
    236368370686157,
    719320896319517,
    393088058075029,
    782048751021037,
    973121922575609
]

medium_small_primes = [
    50020877, 81516553, 15680297, 90581851, 68736721,
    78352691, 50223493, 32618167, 77783243, 75707747,
    85768871, 37239481, 38265803, 88877059, 46509109,
    62098117, 85498649, 25849711, 82554473, 39852167,
    84248993, 96947233, 50138317, 86497559, 73754971,
    79561519, 74793893, 86109053, 98082203, 66439397
]


medium_primes = [
    4021537007,
    8553611153,
    8517840167,
    5784524281,
    1622496451,
    7981809079,
    5771336689,
    2996000281,
    2915992039,
    9562655747,
    2442977083,
    7971934651,
    6824016181,
    1728845597,
    3942616243,
    1716075083,
    6175996831,
    9007594961,
    9660483383,
    2233056877,
    1028802503,
    8289279839,
    7799328241,
    9361518587,
    1809166613,
    2317138123,
    2609085737,
    9366552463,
    5634611303,
    3998121103
]






test_it()