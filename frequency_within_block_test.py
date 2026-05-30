

from __future__ import print_function

import math
from fractions import Fraction
from gamma_functions import *

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

def frequency_within_block_test(bits, alpha=0.01, verbose=False, block_size=None, num_blocks=None):
    # 把序列分成几块，每块看一下 1 的比例是不是偏得离谱
    n = len(bits)
    if block_size is not None:
        M = int(block_size)
        if M <= 0:
            return False, None, {"reason": "块大小得是正数"}
        N = int(math.floor(n/M))
    elif num_blocks is not None:
        N = int(num_blocks)
        if N <= 0:
            return False, None, {"reason": "块数得是正数"}
        M = int(math.floor(n/N))
    else:
        M = 20
        N = int(math.floor(n/M))
        if N > 99:
            N = 99
            M = int(math.floor(n/N))
    
    if len(bits) < 100:
        return False, None, {"reason": "长度不够：这个测试至少要 100 位"}
    if N <= 0 or M <= 0:
        return False, None, {"reason": "分块参数不合适"}
    
    if verbose:
        print("  n = %d" % len(bits))
        print("  N = %d" % N)
        print("  M = %d" % M)
    
    num_of_blocks = N
    block_size = M
    if num_of_blocks * block_size > n:
        return False, None, {"reason": "块数和块大小加起来超过了序列长度"}
    
    proportions = list()
    for i in range(num_of_blocks):
        block = bits[i*(block_size):((i+1)*(block_size))]
        zeroes,ones = count_ones_zeroes(block)
        proportions.append(Fraction(ones,block_size))

    chisq = 0.0
    for prop in proportions:
        chisq += 4.0*block_size*((prop - Fraction(1,2))**2)
    
    p = gammaincc((num_of_blocks/2.0),float(chisq)/2.0)
    success = (p >= alpha)
    return (success, p, None)
