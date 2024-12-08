import cmath
from mpmath import mp

def calc_zeta_by_eta(s, N=1000):
    def eta(s):
        r = 0
        for n in range(1, N, 1):
            r += 1. * (-1)**(n+1) / n **s
        return r
    return eta(s) / (1.-2**(1-s))

mp.dps = 10
for s in [0.5 + 14.134725141*1j, 0.5+1j, 2, 2.00001, 2.000001,  1+1j, 2+2j, 20+20j]:
    #x = get_zeta_val_by_int(s)
    x = calc_zeta_by_eta(s, 1000000)
    print ('calc', s, x)

    x1 = mp.zeta(s)
    print ('should_be', s, x1)
    print ('__')
