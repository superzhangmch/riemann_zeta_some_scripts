from sympy import mobius
from sympy import isprime
from prime_count_directly import count_primes

def is_prime(x):
    return isprime(x)

def Pi(x, adjust=True):
    ret = count_primes(int(x))
    if adjust and type(x) == type(1) and is_prime(x):
        return ret + 0.5
    return ret

def J(x, adjust=True):
    s = 0
    for n in range(1, 1000000000, 1):
        xn = x**(1/n)

        if adjust and type(x) == type(1):
            xn1 = int(round(xn))
            if xn1 **n == x:
                xn = int(xn)

        c = Pi(xn, adjust)
        if c <= 0: break
        s += 1. / n * c
    return s 


def validate_pi_x_mobius(x, adjust=True):
    real = count_primes(x)
    s = 0
    for i in range(1, 1000):
        cur_x_n = x**(1./i) 
        if cur_x_n <= 1.99: break

        if adjust and type(x) == type(1):
            xn1 = int(round(cur_x_n))
            if xn1**i == x:
                cur_x_n = int(cur_x_n)
        s1 = 1. * mobius(i) / i * J(cur_x_n, adjust)
        s += s1
        print (i, round(s1, 3), round(s, 3))
    print (s, 'diff', round(s - real, 3), 'real', real)

if __name__ == '__main__':
    x = 19*19; validate_pi_x_mobius(x, False)
    x = 19; validate_pi_x_mobius(x, True)
    x = 113; validate_pi_x_mobius(x, True)
    x = 123456; validate_pi_x_mobius(x, False)
