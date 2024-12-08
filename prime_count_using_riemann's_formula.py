import mpmath
import numpy as np
from sympy import mobius
import scipy.special as sp
from scipy.special import zeta
import math
import cmath

from scipy.integrate import quad

# https://en.wikipedia.org/wiki/Prime-counting_function

zeta0 = [float(L) for L in open("zeros_100k.txt") if L] # 提前准备的每行一个 zeta 零点。 https://github.com/Sielinski/Riemann-hypothesis/blob/master/zeros_100k
def zt(i):
    return zeta0[i-1]

def f_int_log2(x):
    '''
    int_x^infty (1/[t(t^2-1)ln t]) dt - ln 2
    '''
    def integrand(t): return 1 / (t * (t**2 - 1) * np.log(t))
    end = 100000
    if end < x: end = x * 10
    result, error = quad(integrand, x, end)
    result = -math.log(2.)+result
    return result

#mpmath.mp.dps = 50
#li = mpmath.li
#ei = mpmath.ei
ei = sp.expi
def li(x): return ei(cmath.log(x))
factorial = sp.factorial
#zeta = mpmath.zeta

def pi_x_by_riemann_full(x, N=20, K=500, form=0):
    '''
    N: 莫比乌斯展开项数
    K: 使用的 zeta 零点数
    '''
    result = 0
    for n in range(1, N, 1):
        main_li = li(x**(1./n))
        zero_li = 0
        for k in range(1, K, 1):
            rho0 = 0.5 + 1j * zt(k)
            rho1 = 0.5 - 1j * zt(k)
            zero_li += ei(math.log(x) * rho0 / n)  # 按一般认知：Ei(log(x) * (rho/n)) == Ei(log(x**(rho/n))), 而 Ei(log(x)) == Li(x), 所以按说用 Li(x ** (rho0 / n)) 也行。
            zero_li += ei(math.log(x) * rho1 / n)  # 但是 log(x) 是多值函数。log(x**(rho/n))!= log(x)*rho/n, 所以这里不能简单用 Li(x ** (rho0 / n))而要选对取值分支
                                                   # https://math.stackexchange.com/questions/1144932/curve-profile-for-the-logarithm-integral-sum-term-of-riemann-explicit-formula
                                                   # https://math.stackexchange.com/questions/2066612/evaluating-prime-counting-function-using-riemanns-r-function-and-zeros-of-zeta
        int_log2 = f_int_log2(x) if (form == 0) else 0
        result += 1. * mobius(n) / n * (main_li - zero_li + int_log2)
    if form == 1:
        v = -1 / math.log(x) + 1/math.pi * math.atan(math.pi/math.log(x)) #-1/log(x) + 1/pi * arctan(pi/log(x))
        result += v
    return result

def pi_x_by_mobius_li(x, N=20, form=0):
    # N: 莫比乌斯展开项数
    if form == 0:
        result = 0
        for n in range(1, N, 1):
            result +=  1. * mobius(n) / n * li(x**(1./n))
        return result
    else:
        result = 1
        for n in range(1, N, 1):
            try:
                result += (log(x)**n)/(n*factorial(n)*zeta(n+1))
            except:
                break
        return result

def pi_x_by_li(x):
    # li(x) = int(dt/ln(t)), t=0..x
    return li(x)

def pi_x_by_x_log_x(x):
    # x / log(x)
    return x / math.log(x)

# ============ 
arr = []
for dd in range(1, 20, 1):
    x = 10**dd
    arr1 = ["10e%d" % dd]
    s = pi_x_by_riemann_full(x, N=10, K=500, form=0); arr1.append(s)
    s = pi_x_by_riemann_full(x, N=10, K=500, form=1); arr1.append(s)
    s = pi_x_by_mobius_li(x, N=20, form=0); arr1.append(s)
    s = pi_x_by_mobius_li(x, N=1000, form=1); arr1.append(s)
    s = pi_x_by_li(x);  arr1.append(s)
    s = pi_x_by_x_log_x(x); arr1.append(s)
    arr.append(arr1)
for arr1 in arr:
    print ("|".join(map(str, arr1)))

dd = 20
x = 10**dd
s = pi_x_by_li(x); show_predict(s, dd, m_cnt)
s = pi_x_by_mobius_li(x, N=20, form=0); show_predict(s, dd, m_cnt)
s = pi_x_by_mobius_li(x, N=1000, form=1); show_predict(s, dd, m_cnt)
s = pi_x_by_x_log_x(x); show_predict(s, dd, m_cnt)
s = pi_x_by_riemann_full(x,  N=10, K=500, form=0); show_predict(s, dd, m_cnt)
