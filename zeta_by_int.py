import cmath
from mpmath import mp
import numpy as np
from scipy.special import gamma

def do_path_int(f, arr_path, div_num=1000):
    '''
    离散化计算下沿路径的复积分
    '''
    if type(arr_path) == type([]): arr_path = np.array(arr_path)
    
    integral = 0.0

    # Loop over each consecutive pair of points to form line segments
    for i in range(len(arr_path) - 1):
        z_start = arr_path[i]
        z_end = arr_path[i + 1]
        
        # Calculate the small change in z for each subdivision
        if abs(z_end.real-z_start.real) > 10 or  abs(z_end.imag-z_start.imag) > 10:
            num_subdivisions = 100000
        else:
            num_subdivisions = 1000  # number of subdivisions per line segment
        dz = (z_end - z_start) / num_subdivisions
        
        # Sum up contributions from each sub-segment
        last_v = None
        for j in range(num_subdivisions):
            z0 = z_start + j * dz
            z1 = z_start + (j + 1) * dz
            
            # Approximate integral over the small segment by trapezoidal rule
            v = 0
            try:
                f1, f2 = f(z0), f(z1)
                if f1 is None or f2 is None:
                    v = last_v
                    if v is None:
                        print ('aaaavvvv', i, j, z0, z1, f1, f2)
                else:
                    v = 0.5 * (f1 + f2) * dz
                    last_v = v
            except:
                v = last_v
                print ('bbb', v)
            integral +=v

    return integral

def get_zeta_val_by_int(s):
    def f_int(z):
        try:
            ret = (-z)**s / (cmath.exp(z)-1) / z
            return ret
        except:
            print ('aaaa', z)
            return None
    x = do_path_int(f_int, [500.+1j, -.1 + 1j, -0.1-1j, 500.-1j]) # 用 500 当做 infity
    # 路径积分。从实轴下方的500+1j到上方500-1j之间不需要路径积分 
    zeta_s = gamma(1-s) / (2*cmath.pi*1j) * x 
    return zeta_s

mp.dps = 10
for s in [0.5 + 14.134725141*1j, 2, 2.00001, 2.000001,  1+1j, -2-2j, 2+2j, -20+20j]:
    x = get_zeta_val_by_int(s)
    print ('calc', s, x)

    x1 = mp.zeta(s)
    print ('should_be', s, x1)
    print ('__')
