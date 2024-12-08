from mpmath import mp
import numpy as np
from scipy.special import bernoulli
from scipy.special import loggamma as gammaln

def compute_expression(s, i):
    '''
    直接算阶乘/gamma 函数，要过大溢出。用 log 形式
    '''
    i = float(i)
    
    # Calculate the logarithm of the factorial terms
    log_numerator = gammaln(s + 2*i - 1)  # (s + 2*i - 2)!
    log_denominator = gammaln(2*i + 1) + gammaln(s)  # (2*i)! * (s-1)!
    
    # Compute the expression using the exponential of the difference
    result = np.exp(log_numerator - log_denominator)
    
    return result

def zeta_approximation(s, N=100):
    # Initial sum of the series
    sum_terms = sum(1.0 / n**s for n in range(1, N+1))
    
    # Calculate the correction term using Euler-Maclaurin formula
    correction = 1. / N**(s-1) / (s-1) - 0.5 * N**(-s)
    
    # Bernoulli numbers for the correction terms
    #B = bernoulli(2 * np.arange(1, 10)) # python3
    m = 50
    B = bernoulli(2*m) # bernoulli 数列奇数项取值为0
    
    # Add the Bernoulli correction terms
    correction1 = 0.
    for k in range(1, m+1):
        xxx = compute_expression(s, k)

        dd = 1. * B[2*k] / N**(s-1+2*k) * xxx
        if abs(dd.real) < 1e-20 and abs(dd.imag) < 1e-20: break # 动态决定为达一定精度，要用多少项
        print ('xx', k)
        correction1 += dd
    return sum_terms + correction + correction1

def calc(s):
    mp.dps = 10
    x = mp.zeta(s)
    for N in [100, 1000]:
        approx_zeta = zeta_approximation(s, N=N)
        print ('appr', approx_zeta.real, "\t", approx_zeta.imag)
    print ('should_be', float(x.real), "\t", float(x.imag))

calc(0.5 + 14.134725141*1j)
