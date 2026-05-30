

from __future__ import print_function

import math
import copy
import gf2matrix

def binary_matrix_rank_test(bits, M=32, Q=32, alpha=0.01, verbose=False):
    n = len(bits)
    N = int(math.floor(n/(M*Q))) # 块数
    if verbose:
        print("  块数 %d" % N)
        print("  用掉的 bit 数: %d" % (N*M*Q))
        print("  丢掉的 bit 数: %d" % (n-(N*M*Q))) 
    
    if N < 38:
        return False, None, {"reason": "块数太少：至少要 38 个块才能做这个测试"}
        
    # 下面算的是标准里给的理论概率
    r = M
    product = 1.0
    for i in range(r):
        upper1 = (1.0 - (2.0**(i-Q)))
        upper2 = (1.0 - (2.0**(i-M)))
        lower = 1-(2.0**(i-r))
        product = product * ((upper1*upper2)/lower)
    FR_prob = product * (2.0**((r*(Q+M-r)) - (M*Q)))
    
    r = M-1
    product = 1.0
    for i in range(r):
        upper1 = (1.0 - (2.0**(i-Q)))
        upper2 = (1.0 - (2.0**(i-M)))
        lower = 1-(2.0**(i-r))
        product = product * ((upper1*upper2)/lower)
    FRM1_prob = product * (2.0**((r*(Q+M-r)) - (M*Q)))
    
    LR_prob = 1.0 - (FR_prob + FRM1_prob)
    
    FM = 0      # 满秩矩阵个数
    FMM = 0     # 秩少 1 的矩阵个数
    remainder = 0
    for blknum in range(N):
        block = bits[blknum*(M*Q):(blknum+1)*(M*Q)]
        # 转成矩阵
        matrix = gf2matrix.matrix_from_bits(M,Q,block,blknum) 
        # 算秩
        rank = gf2matrix.rank(M,Q,matrix,blknum)

        if rank == M:
            FM += 1
        elif rank == M-1:
            FMM += 1  
        else:
            remainder += 1

    chisq =  (((FM-(FR_prob*N))**2)/(FR_prob*N))
    chisq += (((FMM-(FRM1_prob*N))**2)/(FRM1_prob*N))
    chisq += (((remainder-(LR_prob*N))**2)/(LR_prob*N))
    p = math.e **(-chisq/2.0)
    success = (p >= alpha)
    
    if verbose:
        print("  满秩矩阵个数  = ", FM)
        print("  秩少 1 的个数 = ", FMM)
        print("  其它情况个数  = ", remainder) 
        print("  卡方值 = ", chisq)

    return (success, p, None)
