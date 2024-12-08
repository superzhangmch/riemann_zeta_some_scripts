import mpmath
import sys
import numpy as np
from sympy import mobius
import scipy.special as sp
from scipy.special import zeta
import math
import cmath

from scipy.integrate import quad
from prime_count_using_mobius import J

# https://en.wikipedia.org/wiki/Prime-counting_function

zeta0 = [float(L.strip()) for L in open("zeros6") if L.strip()] # 提前准备的每行一个 zeta 零点
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

mpmath.mp.dps = 30
li = mpmath.li
ei = mpmath.ei
ei = sp.expi
def li_1(x): return ei(cmath.log(x))

def pi_x_by_riemann_full(x, N=2000, K=1000, form=0):
    '''
    N: 莫比乌斯展开项数
    K: 使用的 zeta 零点数
    '''
    result = 0
    for n in range(1, N, 1):
        cur_x_n = x**(1./n)
        if cur_x_n <= 1.99: break # 只需截止到这里
        if mobius(n) == 0: continue

        use_simple = 0
        if form == 1 and cur_x_n <= 100000000:
            Jx = 1. * mobius(n) / n * J(cur_x_n)
            use_simple = 1
        else:
            main_li = li(x**(1./n))
            zero_li = 0
            for k in range(1, K, 1):
                rho0 = 0.5 + 1j * zt(k)
                rho1 = 0.5 - 1j * zt(k)
                zero_li += ei(math.log(x) * rho0 / n)
                zero_li += ei(math.log(x) * rho1 / n)
            int_log2 = f_int_log2(x)
            Jx = 1. * mobius(n) / n * (main_li - zero_li + int_log2)
        result += Jx
        #print ('nn', n, Jx, result, use_simple, int(cur_x_n), type(result))
    return result

if __name__ == '__main__':
    from cnt import m_cnt
    if 1: # 看看黎曼素数计数公式准确度如何. 可以看出当x是素数的时候，确实是累加了0.5，x+epsilon后才变为+1
        for i in range(1, 2000):
            x = 0.01 * i
            s = pi_x_by_riemann_full(x,  N=500, K=1000, form=0)
            print ('final', x, s)
    if 0:
        for i in range(1, 20):
            dd = i
            x = 10**dd
            s0 = pi_x_by_riemann_full(x,  N=500, K=1000, form=0) # 所有项都zeta零点参与算
            s1 = pi_x_by_riemann_full(x,  N=500, K=1000, form=1) # 只头几项用zeta零点参与算
            real = m_cnt[dd]
            print (dd, s0, s1, abs(s0-real), abs(s1-real))
