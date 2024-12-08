from mpmath import mp # 他可以高精度计算 zeta(s) 的值

def zeta_approximation(s, N=100):
    sum_terms = sum(1.0 / n**s for n in range(1, N+1))
    return sum_terms

mp.dps = 10
s = -0.1 + 49.773832478 * 1j  # Example complex number
s = 1 + 1.0 * 1j  # Example complex number
x = mp.zeta(s)
approx_zeta = zeta_approximation(s, N=1000000)
print 'appr', approx_zeta.real, "\t", approx_zeta.imag
print 'real', float(x.real), "\t", float(x.imag)

s = 0.5 + 1.0 * 1j
for N in [100, 1000, 10000, 100000, 1000000]:
    approx_zeta = zeta_approximation(s, N=N)
    print 'appr', approx_zeta.real, "\t", approx_zeta.imag

# 有解析延拓的 zeta(s) 的真实取值
mp.dps = 10
x = mp.zeta(s)
print 'should_be', float(x.real), "\t", float(x.imag)
