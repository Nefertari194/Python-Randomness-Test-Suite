

from __future__ import print_function

import math

# 这里就是简单数一下 0/1 的个数
def count_ones_zeroes(bits):
    ones = 0
    zeroes = 0
    for bit in bits:
        if (bit == 1):
            ones += 1
        else:
            zeroes += 1
    return (zeroes, ones)

def runs_test(bits, alpha=0.01, verbose=False):
    n = len(bits)
    if n == 0:
        return False, None, {"reason": "数据是空的"}

    zeroes, ones = count_ones_zeroes(bits)

    prop = float(ones)/float(n)
    if verbose:
        print("  prop ", prop)

    tau = 2.0/math.sqrt(n)
    if verbose:
        print("  tau ", tau)

    if abs(prop-0.5) > tau:
        return False, None, {"reason": "0/1 比例偏得有点多，这个测试按标准就不做了"}

    vobs = 1.0
    for i in range(n-1):
        if bits[i] != bits[i+1]:
            vobs += 1.0

    if verbose:
        print("  vobs ", vobs)
      
    p = math.erfc(abs(vobs - (2.0*n*prop*(1.0-prop)))/(2.0*math.sqrt(2.0*n)*prop*(1-prop) ))
    success = (p >= alpha)
    return (success, p, None)
