

from __future__ import print_function

import math

def count_ones_zeroes(bits):
    ones = 0
    zeroes = 0
    for bit in bits:
        if (bit == 1):
            ones += 1
        else:
            zeroes += 1
    return (zeroes, ones)

def monobit_test(bits, alpha=0.01, verbose=False):
    n = len(bits)
    if n == 0:
        return False, None, {"reason": "数据是空的"}
    
    zeroes, ones = count_ones_zeroes(bits)
    s = abs(ones-zeroes)
    if verbose:
        print("  1 的个数 = %d" % ones)
        print("  0 的个数 = %d" % zeroes)
    
    p = math.erfc(float(s)/(math.sqrt(float(n)) * math.sqrt(2.0)))
    
    success = (p >= alpha)
    return (success, p, None)
    
