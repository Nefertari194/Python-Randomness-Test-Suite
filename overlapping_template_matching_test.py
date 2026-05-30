

from __future__ import print_function

import math
from gamma_functions import *

def lgamma(x):
    return math.log(gamma(x))
    
def Pr(u, eta):
    if ( u == 0 ):
        p = math.exp(-eta)
    else:
        sum = 0.0
        for l in range(1,u+1):
            sum += math.exp(-eta-u*math.log(2)+l*math.log(eta)-lgamma(l+1)+lgamma(u)-lgamma(l)-lgamma(u-l+1))
        p = sum
    return p

def overlapping_template_matching_test(bits, blen=6, alpha=0.01, verbose=False, block_size=None, num_blocks=None, allow_degraded=False):
    n = len(bits)
    
    m = int(blen)
    if m < 2:
        m = 2
    # 模板（这里按原思路用全 1）
    B = [1 for x in range(m)]
    
    K = 5
    M = 1062
    N = 968
    if allow_degraded:
        if block_size is not None:
            M = int(block_size)
            if M <= 0:
                return False, None, {"reason": "块大小得是正数"}
            N = int(math.floor(n / M))
        elif num_blocks is not None:
            N = int(num_blocks)
            if N <= 0:
                return False, None, {"reason": "块数得是正数"}
            M = int(math.floor(n / N))
        else:
            N = max(1, int(math.floor(n / 128)))
            M = int(math.floor(n / N))

        if N <= 0 or M <= 0:
            return False, None, {"reason": "分块参数不合适"}
        if M <= m:
            return False, None, {"reason": "块太短了，模板放不下"}
        if N < 2:
            return False, None, {"reason": "块数太少，结果不太靠谱"}
    else:
        if len(bits) < (M*N):
            return False, None, {"reason": "长度不够：至少要 %d 位（现在 %d）" % (M*N, len(bits))}
    
    blocks = list() # 分成若干块
    for i in range(N):
        blocks.append(bits[i*M:(i+1)*M])

    # 统计每块里命中次数的分布
    v=[0 for x in range(K+1)] 
    for block in blocks:
        count = 0
        for position in range(M-m+1):
            if block[position:position+m] == B:
                count += 1
            
        if count >= (K):
            v[K] += 1
        else:
            v[count] += 1

    chisq = 0.0  # 算卡方
    pi = [0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865] # 参考参数
    piqty = [int(x*N) for x in pi]
    
    lambd = (M-m+1.0)/(2.0**m)
    eta = lambd/2.0
    sum = 0.0
    for i in range(K): # 算概率
        pi[i] = Pr(i, eta)
        sum += pi[i]

    pi[K] = 1 - sum;

    sum = 0    
    chisq = 0.0
    for i in range(K+1):
        chisq += ((v[i] - (N*pi[i]))**2)/(N*pi[i])
        sum += v[i]
        
    p = gammaincc(5.0/2.0, chisq/2.0) # 算 p 值

    if verbose:
        print("  B = ",B)
        print("  m = ",m)
        print("  M = ",M)
        print("  N = ",N)
        print("  K = ",K)
        print("  model = ",piqty)
        print("  v[j] =  ",v) 
        print("  chisq = ",chisq)
    
    success = (p >= alpha)
    return (success, p, None)
