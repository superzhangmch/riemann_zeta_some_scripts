# 从 https://github.com/RobinLinus/prime-counting-function/blob/master/prime-counting.js 逐行翻译到了python

def sieve_of_eratosthenes(limit):
    # 初始化一个布尔数组，其值初始化为True
    prime = [True for _ in range(limit+1)]
    p = 2
    while (p * p <= limit):
        # 如果prime[p]没有被改变，那么它一定是一个素数
        if (prime[p] == True):
            # 更新所有p的倍数为非素数
            for i in range(p * p, limit+1, p):
                prime[i] = False
        p += 1
    # 收集所有素数
    prime_numbers = []
    for p in range(2, limit):
        if prime[p]:
            prime_numbers.append(p)
    return prime_numbers

Primes = sieve_of_eratosthenes(1300000)
MAX_PRIME = Primes[-1]  # The largest prime in our database
MAX_N = MAX_PRIME * MAX_PRIME  # Our database enables us to count up to this number
#print (MAX_PRIME, MAX_N)

def count_primes(n):
    if n < 2:
        return 0
    if n == 2:
        return 1
    if n == 3:
        return 2
    if n > MAX_N:
        raise ValueError('N is too big. Bigger database of primes required!')
    
    count_primes = n  # We assume all numbers are prime and then we subtract the composites
    index = 0
    p = Primes[0]

    while p * p <= n:  # We iterate over all primes up to sqrt(n)
        m = n // p  # FIXME: this integer division could be done faster
        count_primes -= count_almost_primes_cached(m, p, index)  # subtract all composites which's smallest factor is p

        index += 1
        p = Primes[index]
    
    return count_primes

cache2 = {}

def count_almost_primes_cached(m, p, index):
    # We know these results without computation:
    if p == 2:
        return m  # All numbers are 2-almost-primes
    if m == 0:
        return -1  # FIXME: where are these artifacts coming from??
    if m < p:
        return 0  # No number has a factor bigger than itself
    if m == p or m - 1 == p:
        return 1  # There is only one number with a factor as big as itself (or one less)
    
    if p * p > m:  # This means, all "p-almost-primes" up to m are actually primes
        # So let's count only the primes
        
        if m > MAX_PRIME:
            return count_primes(m) - index  # m is still too large. Let's recurse
        
        # Yey, we can calculate count_primes(m) "by hand"!
        count_primes_m = index
        while Primes[count_primes_m] <= m:
            count_primes_m += 1  # FIXME: Binary search here
        return count_primes_m - index

    # Otherwise, there's actually work to do. Let's cache it
    cache_key = f"{m}_{p}"
    if cache_key not in cache2:
        cache2[cache_key] = _count_almost_primes(m, p)
    return cache2[cache_key]

def _count_almost_primes(m, p):
    index = 0
    q = Primes[0]

    almost_primes = m  # Assume all numbers are p-almost-prime and then subtract the rest
    while q < p:  # We iterate over all primes smaller than p
        m1 = m // q  # FIXME: this integer division could be done faster
        almost_primes -= count_almost_primes_cached(m1, q, index)

        index += 1
        q = Primes[index]
    
    almost_primes -= index
    return almost_primes

#for i in range(1, 9): print (i, count_primes(10**i))
