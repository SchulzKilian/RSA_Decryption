import math


def compute_primes(x):
    length = int(math.sqrt(x)) +1
    l = [False for i in range(length-1)]
    x = True
    for ind,i in enumerate(l):
        if l[ind] or ind==0:
            continue
        numb = ind+1
        result = x/numb
        if result == int(result):
            return result, numb
        else:
            l[ind]= True
            k = numb
            while k*numb <= length-1 and k<10:
                l[k*numb-1]= True
                k += numb


print(compute_primes(22742734291*52711))