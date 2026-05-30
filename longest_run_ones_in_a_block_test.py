

from __future__ import print_function

import math
from gamma_functions import *

def probs(K,M,i):
    M8 =      [0.2148, 0.3672, 0.2305, 0.1875]
    M128 =    [0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124]
    M512 =    [0.1170, 0.2460, 0.2523, 0.1755, 0.1027, 0.1124]
    M1000 =   [0.1307, 0.2437, 0.2452, 0.1714, 0.1002, 0.1088]
    M10000 =  [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]
    if (M == 8):        return M8[i]
    elif (M == 128):    return M128[i]
    elif (M == 512):    return M512[i]
    elif (M == 1000):   return M1000[i]
    else:               return M10000[i]

def longest_run_ones_in_a_block_test(bits, alpha=0.01, verbose=False):
    n = len(bits)

    if n < 128:
        return False, None, {"reason": "长度不够：这个测试至少要 128 位"}
    elif n<6272:
        M = 8
    elif n<750000:
        M = 128
    else:
        M = 10000
            
    # 根据块大小确定几个参数
    if M==8:
        K=3
        min_N = 16
    elif M==128:
        K=5
        min_N = 49
    else:
        K=6
        min_N = 75

    N = int(math.floor(n/M))
    if N < min_N:
        return False, None, {"reason": "有效块数太少（N=%d < %d）" % (N, min_N)}
        
    # 统计表
    v = [0,0,0,0,0,0,0]

    for i in range(N): # 遍历每个块
        # 找最长连续 1
        block = bits[i*M:((i+1)*M)]
        
        run = 0
        longest = 0
        for j in range(M): # 数位
            if block[j] == 1:
                run += 1
                if run > longest:
                    longest = run
            else:
                run = 0

        if M == 8:
            if longest <= 1:    v[0] += 1
            elif longest == 2:  v[1] += 1
            elif longest == 3:  v[2] += 1
            else:               v[3] += 1
        elif M == 128:
            if longest <= 4:    v[0] += 1
            elif longest == 5:  v[1] += 1
            elif longest == 6:  v[2] += 1
            elif longest == 7:  v[3] += 1
            elif longest == 8:  v[4] += 1
            else:               v[5] += 1
        else:
            if longest <= 10:   v[0] += 1
            elif longest == 11: v[1] += 1
            elif longest == 12: v[2] += 1
            elif longest == 13: v[3] += 1
            elif longest == 14: v[4] += 1
            elif longest == 15: v[5] += 1
            else:               v[6] += 1
    
    # 算卡方
    chi_sq = 0.0
    for i in range(K+1):
        p_i = probs(K,M,i)
        upper = (v[i] - N*p_i)**2
        lower = N*p_i
        chi_sq += upper/lower
    if verbose:
        print("  n = "+str(n))
        print("  K = "+str(K))
        print("  M = "+str(M))
        print("  N = "+str(N))
        print("  chi_sq = "+str(chi_sq))
    p = gammaincc(K/2.0, chi_sq/2.0)
    
    success = (p >= alpha)
    return (success, p, None)
