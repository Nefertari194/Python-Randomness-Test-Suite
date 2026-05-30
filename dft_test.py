

from __future__ import print_function

import math
import numpy
import sys

def dft_test(bits, alpha=0.01, verbose=False):
    n = len(bits)
    if (n % 2) == 1:        # 这里必须是偶数长度，不够的话就丢掉最后 1 位
        bits = bits[:-1]
        n = len(bits)

    ts = list()             # 把 0/1 变成 -1/+1
    for bit in bits:
        ts.append((bit*2)-1)

    ts_np = numpy.array(ts)
    fs = numpy.fft.fft(ts_np)  # 做离散傅里叶变换
   
    if sys.version_info > (3,0):
        mags = abs(fs)[:n//2] # 只用前半段的幅值
    else:
        mags = abs(fs)[:n/2] # 只用前半段的幅值
    
    T = math.sqrt(math.log(1.0/0.05)*n) # 上阈值
    N0 = 0.95*n/2.0
    if verbose:
        print("  N0 = %f" % N0)

    N1 = 0.0   # 统计小于阈值的点（按原实现）
    for mag in mags:
        if mag < T:
            N1 += 1.0
    if verbose:
        print("  N1 = %f" % N1)
    d = (N1 - N0)/math.sqrt((n*0.95*0.05)/4) # 算 p 值
    p = math.erfc(abs(d)/math.sqrt(2))

    success = (p >= alpha)
    return (success, p, None)
