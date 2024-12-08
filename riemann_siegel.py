from sympy import symbols, diff
import sympy
import math
import cmath
import copy
import numpy as np
import mpmath
from scipy.special import gamma, loggamma
from scipy.special import bernoulli

# <<
def get_theta_coefficients():
    # https://arxiv.org/pdf/2004.00926
    B = bernoulli(2*50)
    coef = []
    for k in range(1, 50):
        aa = (2**(2*k-1)-1) * B[2*k]
        bb = 4*k*(2*k-1) * 2**(2*k-1)
        cc = aa/bb
        coef.append(cc) 
    return coef

def theta_t(t):
    '''
    https://en.wikipedia.org/wiki/Riemann%E2%80%93Siegel_theta_function
    https://mathworld.wolfram.com/Riemann-SiegelFunctions.html
    https://arxiv.org/pdf/2004.00926
    = arg(gamma(1/4+1/2*t*1j))- 1/2*t*log(pi)
    '''
    log = math.log
    pi = math.pi
    # 还有个 arctan(e^{−πt}) 项，在 t 稍大后，可以省去
    # 下面 /t, /t**3, /t**5系数可用get_theta_coefficients() 获得
    return -t/2*log(2*pi/t) - t/2. - pi/8 + 1/48./t + 7/5760./t**3 + 31/80640./t**5


def theta_t_1(t):
    '''
    按说用它就行了，但是需要返回的gamma的辐角能保持连续, 下面这样做不到。故简单这样做不行
    '''
    return cmath.phase(gamma(1./4+t*1/2j)) - 0.5*t*math.log(math.pi)

if 0:
    for i in range(1000, 5000):
        print (i, theta_t_1(0.1*i), theta_t(0.1*i))
    t = 10000
    print (theta_t(t), theta_t_1(t))
# >>

def get_psi_func(rt_cnt=5):
    # https://mathworld.wolfram.com/Riemann-SiegelFormula.html

    assert rt_cnt <= 5
    P = symbols('P')

    cos = sympy.cos
    pi = sympy.pi

    psi_expr = cos(2*pi*(P*P-P-sympy.Rational(1, 16))) / cos(2*pi*P)
    
    drctv_n = {5:16, 4:13, 3: 10, 2: 7, 1: 4, 0: 1}[rt_cnt]
    psi_d_arr = [psi_expr]
    for n in range(1, drctv_n, 1):
        print ('calc directive', n)
        d = diff(psi_expr, P, n)
        psi_d_arr.append(d)

    def get_f(d):
        def f(p_v):
            v = d.subs(P, p_v).evalf()
            return v
        return f

    psi_d = []
    for d in psi_d_arr:
        psi_d.append(get_f(d))

    pi = math.pi
    psi_d[0] = lambda p: math.cos(2*pi*(p*p-p-1/16.)) / math.cos(2*pi*p)

    return psi_d

psi_d = []
def Rt(t, rt_cnt=5):
    pi = math.pi
    t2pi = t / (2*pi)
    N = int(math.sqrt(t2pi))
    p = math.sqrt(t2pi)-N

    rt = (-1)**(N-1) * t2pi**(-1/4.)
    s = 0
    for k in range(rt_cnt):
        if k == 0:
            ck = psi_d[0](p)
        elif k == 1:
            ck = -psi_d[3](p)/(96*pi**2)
        elif k == 2:
            ck = psi_d[2](p)/(64*pi**2)+psi_d[6](p)/(18432*pi**4)
        elif k == 3:
            ck = -psi_d[1](p)/(64*pi**2)-psi_d[5](p)/(3840*pi**4)-psi_d[9](p)/(5308416*pi**6) 
        elif k == 4:
            ck = psi_d[0](p)/(128*pi**2)+19*psi_d[4](p)/(24576*pi**4)+11*psi_d[8](p)/(5898240*pi**6)+psi_d[12](p)/(2038431744*pi**8)
        elif k == 5:
            ck = -5*psi_d[3](p)/(3072*pi**4)-901*psi_d[7](p)/(82575360*pi**6)-7*psi_d[11](p)/(849346560*pi**8)-psi_d[15](p)/(978447237120*pi**(10))
        else:
            assert False
        ck1 = ck * t2pi**(-k/2.)
        #print ('k', k, ck, ck1)
        s += ck1
    return rt * s

def Z_t(t, rt_cnt=5):
    # https://mathworld.wolfram.com/Riemann-SiegelFormula.html
    # https://mathworld.wolfram.com/Riemann-SiegelFunctions.html
    theta_v = theta_t(t)
    t2pi = t / (2*math.pi)
    vt = int(math.sqrt(t2pi))
    #print ('loop cnt', vt)
    s = 0
    for k in range(1, vt+1, 1):
        s1 = 1 / k**0.5 * math.cos(theta_v - t * math.log(k))
        #if k % 100000 == 0: print ('kk', k, s1, s)
        s += s1
    s *= 2

    rt = Rt(t, rt_cnt)
    z_v = s + rt

    # Z(t) = exp(I*theta(t)) * zeta(1/2+t*I)
    zeta_v =  z_v / cmath.exp(theta_v * 1j)
    return s + rt, zeta_v

def check(t):
    print ('t', t, 'result', Z_t(t))
    print ('real value', mpmath.zeta(0.5+1j * t))
    print ('--')

rt_cnt = 5
psi_d = get_psi_func(rt_cnt)

check(14.134725142) # zero point
check(74920.827498994)  # zero point
check(1132490.658714411)  # zero point
check(10000000001.0405584463652298444148667285526) #  # zero point
check(10e10)
for t in [1, 10, 100, 1000, 10000, 10e5, 10e6, 10e7, 10e10, 10e12, 10e15, 10e20]:
    b, c, d, e, f= Rt(t, 1), Rt(t, 2), Rt(t, 3),Rt(t, 4), Rt(t, 5)
    arr = [t, b, c, c-b, d,d-c, e, e-d, f, f-e]
    print ("|", "|".join(map(str, arr)), "|")

# 找几个零点
last_v = 1 
print ('---')
for i in range(1000, 1000 * 5):
    t = i * 0.01
    a, b = Z_t(t, rt_cnt=2)
    if last_v * a < 0: # Z_t(t) 是实值的，看符号变化找零点
        print (t) 
    last_v = a 
    
'''
# 找到的几个零点
10.0 # 跳过这个
14.14
21.03
25.02
30.43
32.94
37.59
40.92
43.33
48.01
49.78
'''
